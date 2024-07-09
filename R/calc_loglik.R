# Function to turn parameters into a string
pars_to_string <- function(pars) {

  # Paste together the parameter names and values
  paste(paste(names(pars), as.list(pars), sep = " = "), collapse = ", ")

}

# Function to calculate the likelihood of the data given the model
calc_loglik <- function(

  data, pars, island_age, M, nmax,
  control = list(), optimized = FALSE

) {

  # data: list of clade data
  # pars: model parameters
  # island_age: the age of the island
  # M: size of the mainland pool
  # nmax: maximum number of unobserved species allowed per clade
  # control: control parameters for the integrator
  # optimized: is this function called from within an optimizer?

  # TODO: Here the original code tells the user if any parameter is infinite.

  # Turn the parameters back into a list if not already
  pars <- as.list(pars)

  # Check the parameters (they will have been changed by the optimizer)
  tryCatch({ check_par_values(pars) }, error = function(error) {

    # Preamble if the problem is encountered during maximization
    if (optimized) error <- paste("During optimization,", error)

    # Append parameter values
    error <- paste0(error, "\nParameters: ", pars_to_string(pars))

    # Error
    stop(error)

  })

  # Compute the likelihood of a clade leaving no descendants
  loglik <- integrate_clade(island_age, pars, nmax, control)

  # TODO: Correct the likelihood of no extand descendant by a mainland clade
  # by a conditioning factor. This is logcond in the original code, which
  # is computed by the function logcondprob(), taking as arguments the
  # type of conditioning, the number of clades on the island and the
  # probability of no descendants we just computed. The conditioning can be
  # either on island age only, island age and extant island biota, or island
  # age and at least a certain number of colonizations (last option not
  # available for the island-wide model).

  # How many of the mainland pool made it to the island?
  ncol <- length(data)

  # How many did not?
  not_present <- M - ncol

  # Sum it over all mainland species that did not colonize
  loglik <- not_present * loglik

  # For each colonist clade on the island...
  for (i in 1:ncol) {



    # Incorporate its likelihood given the parameters
    loglik <- loglik +
      with(data[[i]], integrate_clade(
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

  # TODO:

  # If the likelihood could not be computed
  # if (is.na(loglik)) {
  #
  #   # Say it
  #   message("Some parameter values caused numerical problems.")
  #
  #   # Equate it to zero
  #   return(-Inf)
  #
  # }

  # TODO: Hard to trigger that one.

  return(loglik)

}
