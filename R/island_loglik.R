# Function to turn parameters into a string
pars_to_string <- function(pars) {

  # Paste together the parameter names and values
  paste(paste(names(pars), as.list(pars), sep = " = "), collapse = ", ")

}

# Function to check the parameters during the likelihood calculation
check_pars_loglik <- function(pars, optimized) {

  # pars: the parameters
  # optimized: is this function called from within an optimizer?

  # Wrapper around the parameter checker with custom error message
  tryCatch({ check_par_values(pars) }, error = function(error) {

    # Preamble if the problem is encountered during maximization
    if (optimized) error <- paste("During optimization,", error)

    # Append parameter values
    error <- paste0(error, "\nParameters: ", pars_to_string(pars))

    # Error
    stop(error)

  })
}

# Function to compute the approximate likelihood of no descendants
approximate_loglik0 <- function(pars, island_age) {

  # pars: model parameters
  # island_age: the island age

  # Compute
  with(pars, -log(mu + gamma) + log(mu + gamma * exp(-(mu + gamma) * island_age)))

}

# Function to condition the likelihood on a certain number of colonizations
condition_loglik <- function(loglik, nobs, ncond, logp0) {

  # loglik: the log-likelihood thus far
  # nobs: the observed number of colonizers
  # ncond: the number of colonizations to condition on
  # logp0: the probability of no descendant for a mainland clade

  # Note: conditioning on zero colonizations means only conditioning on island
  # age, while numbers greater than zero mean conditioning on non-extinction
  # of the island biota and at least a certain number of colonization events
  # (more than one is not yet available for the island-wide model).

  # Safety check
  testit::assert(is_number(loglik))
  testit::assert(is_positive_integer(nobs))
  testit::assert(is_positive_integer(ncond, strict = TRUE))
  testit::assert(is_number(logp0))

  # Also that
  testit::assert(ncond <= nobs)

  # TODO: Maybe an error message here if needed?

  # Create some bounds (factor is arbitrary here)
  nmax <- min(nobs, 2L * ncond)

  # Initalize the amount to substract to the log-likelihood
  logcond <- 0

  # Check whether the probability of no descendant is almost zero
  is_nearly_zero <- exp(logp0) == 1 & logp0 < 0

  # Note: numerically it is possible for a negative number to still have its
  # exponential evaluated as exactly one, if it is very close to zero.

  # Not the probability of no descendants
  lognotp0 <- if (is_nearly_zero) log(-logp0) else log1p(-exp(logp0))

  # Prepare a container for clade probabilities
  logpc <- rep(0, nmax + 1L)

  # Range of indices
  ii <- 0:nmax

  # Compute those clade (log-)probabilities
  logpc[ii + 1] <- lgamma(nobs + 1) - lgamma(ii + 1) - lgamma(nobs - ii + 1) +
    (nobs - ii) * logp0 + ii * lognotp0

  # Convert them back into probabilities
  pc <- exp(logpc)

  # Key indices
  jj <- (ncond + 1):(nmax + 1)

  # Special case of a sum of probabilities smaller than one
  if (sum(pc) < 1) return(loglik - log1p(-sum(pc[-jj])))

  # Otherwise, issue a message
  message("An approximation of the likelihood-conditioning must be made. Results may be unreliable.")

  # TODO: If verbose?

  # And make an approximation
  return(loglik - log(sum(pc[jj])))

}

# Function to calculate the likelihood of the data given the model
island_loglik <- function(

  data, pars, island_age, M, nmax, control = list(), ncond = 0L,
  optimized = FALSE

) {

  # data: list of clade data
  # pars: model parameters
  # island_age: the age of the island
  # M: size of the mainland pool
  # nmax: maximum number of unobserved species allowed per clade
  # control: control parameters for the integrator
  # ncond: the number of colonization events to condition on
  # optimized: is this function being called from within an optimizer?

  # Turn the parameters back into a list if not already
  pars <- as.list(pars)

  # Check the parameters (they will have been changed by the optimizer)
  check_pars_loglik(pars, optimized)

  # TODO: Here the original code tells the user if any parameter is infinite.

  # Compute the likelihood of a clade leaving no descendants
  logp0 <- clade_loglik(island_age, pars, nmax, control)

  # Approximate the likelihood of no descendants if needed
  if (logp0 >= 0 & with(pars, mu / lambda_c > 100))
    logp0 <- with(pars, approximate_logp0(pars, island_age))

  # TODO: This might not trigger if the likelihood function breaks from within,
  # which will happen if a positive likelihood is encountered.

  # How many of the mainland pool made it to the island?
  ncol <- length(data)

  # How many did not?
  not_present <- M - ncol

  # Sum it over all mainland species that did not colonize
  loglik <- not_present * logp0

  # Condition the likelihood if needed
  if (ncond > 0L) loglik <- condition_loglik(loglik, ncol, ncond, logp0)

  # For each colonist clade on the island...
  for (i in 1:ncol) {

    # Incorporate its likelihood given the parameters
    loglik <- loglik +
      with(data[[i]], clade_loglik(
        island_age,
        pars,
        nmax,
        is_present, tcol,
        branching_times, tmax, tmin,
        control
      ))

  }

  # TODO: The likelihood function switches to a different calculation if the
  # rate of cladogenesis becomes infinite, with a message (somehow this
  # is only done if there are zero missing species).

  return(loglik)

}
