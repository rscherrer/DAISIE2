# TODO: Make sure default parameters are only where needed.

# Function to calculate the likelihood of the data given the model
calc_loglik <- function(data, pars, island_age, M, nmax) {

  # data: list of clade data
  # pars: model parameters
  # island_age: the age of the island
  # M: size of the mainland pool
  # nmax: maximum number of unobserved species allowed per clade

  # Turn the parameters back into a list if not already
  pars <- as.list(pars)

  # Compute the likelihood of a clade leaving no descendants
  loglik <- integrate_clade(island_age, pars, nmax)

  # TODO: Don't forget to pass extra arguments too.

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
        branching_times, tmin, tmax
      ))

  }

  return(loglik)

}
