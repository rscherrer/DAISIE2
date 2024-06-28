# Function to return names of model parameters
par_names <- function() {

  # Parameter names
  return(c("lambda_c", "mu", "gamma", "lambda_a"))

}

# Check parameters passed to a function
check_pars <- function(pars) {

  # pars: the parameters

  # Must be a list
  testit::assert(is.list(pars))

  # Parameters must be named
  testit::assert(!is.null(names(pars)))

  # Check that all the expected parameters are there
  testit::assert(all(names(pars) %in% par_names()))
  testit::assert(all(par_names() %in% names(pars)))

  # Check that all parameters are positive numbers
  testit::assert(is_number(unlist(pars), scalar = FALSE, sign = 1))

}

# Function to check that the timings are good
check_times <- function(tcol, tmin, tmax, branching_times, island_age) {

  # tcol: known colonization time (could be NULL)
  # tmin: minimum colonization time (could be NULL)
  # tmax: maximum colonization time (could be NULL)
  # branching_times: vector of branching times (could be empty)
  # island_age: age of the island

  # There must be colonization otherwise there is no clade
  testit::assert(!is.null(tcol) | !is.null(tmax))

  # If known colonization time, upper and lower bounds should not be set
  if (!is.null(tcol)) testit::assert(is.null(tmax) & is.null(tmin))

  # If upper bound given, there should be no known colonization time
  if (!is.null(tmax)) testit::assert(is.null(tcol))

  # If a lower bound is given, then an upper bound should have been given too
  if (!is.null(tmin)) testit::assert(!is.null(tmax))

  # Check that the colonization times are all negative numbers
  if (!is.null(tcol)) testit::assert(is_number(tcol, sign = -1))
  if (!is.null(tmax)) testit::assert(is_number(tmax, sign = -1))
  if (!is.null(tmin)) testit::assert(is_number(tmin, sign = -1))

  # And that if provided branching times are too
  if (length(branching_times) > 0L)
    testit::assert(all(is_number(branching_times, scalar = FALSE, sign = -1)))

  # Make sure that timings are in chronological order
  times <- c(island_age, tmax, tcol, tmin, branching_times)
  testit::assert(all(diff(times) > 0))

}

# Function to check a clade
check_clade <- function(clade, island_age = -Inf) {

  # clade: the clade
  # island_age: age of the island

  # Check that island age is a negative number
  testit::assert(is_number(island_age, sign = -1))

  # Must be a list
  testit::assert(is.list(clade))

  # Must contain five elements
  testit::assert(length(clade) == 5L)

  # Items must have names
  testit::assert(!is.null(names(clade)))

  # Must contain the following items
  clade_items <- c("is_present", "tcol", "branching_times", "tmax", "tmin")
  testit::assert(all(names(clade) %in% clade_items))
  testit::assert(all(clade_items %in% names(clade)))

  # TODO: Should we put item names in a separate function? As for parameters?

  # Make sure whether the mainland colonist is present is a yes or no
  with(clade, {

    # TODO: Maybe create a function for yes/no

    testit::assert(is.logical(clade$is_present))
    testit::assert(length(clade$is_present) == 1L)
    testit::assert(!is.na(clade$is_present))

  })

  # Check the timings
  with(clade, check_times(tcol, tmin, tmax, branching_times, island_age))

}

# Function to check the input data
check_data <- function(data) {

  # data: the list of clades

  # Must be a list
  testit::assert(is.list(data))

  # Must have at least one colonist clade
  testit::assert(length(data) > 0L)

  # TODO: What about an empty island?

  # Each element in the list must be a clade
  for (i in seq(data)) check_clade(data[[i]])

}

# Function to find maximum likelihood estimates
daisie_ml <- function(

  data, pars, island_age, M, nmax,
  optimmethod = "subplex",
  tol = c(1e-4, 1e-5, 1e-7),
  maxiter = 1000L * round((1.25)^length(pars)),
  num_cycles = 1L, jitter = 0

) {

  # data: the input data
  # pars: model parameters
  # island_age: age of the island
  # M: size of the mainland pool
  # nmax: maximum allowed number of unobserved species
  # tol, maxiter, num_cycles, jitter: arguments for the optimizer

  # TODO: Actually it's c(tol, maxiter) that is an argument.

  # TODO: Handle methodological arguments better (pass those for the integrator).

  # Check the input data
  check_data(data)

  # Check the parameters
  check_pars(pars)

  # Check that island age is a negative number (i.e. time ago)
  testit::assert(is_number(island_age, scalar = TRUE, sign = -1))

  # Check that the mainland pool is a strictly positive number
  testit::assert(is_number(M, scalar = TRUE, integer = TRUE, sign = 1))

  # TODO: What if it is zero?

  # Check that the maximum number of unobserved species is a positive number
  testit::assert(is_number(nmax, scalar = TRUE, integer = TRUE, sign = 1))

  # We must turn the parameters into a vector for DDD::optimizer
  pars <- unlist(pars)

  simplex

  # Perform the likelihood search
  optimizer(

    # The likelihood function
    fun = calc_loglik,

    # Vector of initial guesses for parameter values
    trparsopt = pars,

    # Other arguments of the likelihood function
    island_age = island_age,
    M = M,
    nmax = nmax,

    # Parameters of the optimizer
    optimmethod = optimmethod,
    optimpars = c(tol, maxiter),
    num_cycles = num_cycles,
    jitter = jitter

  )
}
