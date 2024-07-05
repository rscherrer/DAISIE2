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
  testit::assert(is_positive_vector(unlist(pars)))

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
  if (!is.null(tcol)) testit::assert(is_negative(tcol))
  if (!is.null(tmax)) testit::assert(is_negative(tmax))
  if (!is.null(tmin)) testit::assert(is_negative(tmin))

  # And that if provided branching times are too
  if (length(branching_times) > 0L)
    testit::assert(all(is_negative_vector(branching_times)))

  # Make sure that timings are in chronological order
  times <- c(island_age, tmax, tcol, tmin, branching_times)
  testit::assert(all(diff(times) > 0))

}

# Function to check a clade
check_clade <- function(clade, island_age = -Inf) {

  # clade: the clade
  # island_age: age of the island

  # Check that island age is a negative number
  testit::assert(is_negative(island_age))

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

  # Make sure whether the mainland colonist is present is a yes or no
  testit::assert(is_yes_no(clade$is_present))

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
  control_ml = list(),
  control_ode = list(),
  method = "subplex",
  ncycles = 1L,
  tol = 1e-6

) {

  # data: the input data
  # pars: model parameters
  # island_age: age of the island
  # M: size of the mainland pool
  # nmax: maximum allowed number of unobserved species
  # control_ml: named list of options for the likelihood maximization algorithm
  # control_ode: named list of options for the integrator
  # method: which method (of "simplex" or "subplex") to use for optimization?
  # ncycles: number of optimization cycles
  # tol: tolerance criterion between cycles

  # Check the input data
  check_data(data)

  # Check the parameters
  check_pars(pars)

  # Check that island age is a negative number (i.e. time ago)
  testit::assert(is_negative(island_age))

  # Check that the mainland pool is a strictly positive number
  testit::assert(is_positive_integer(M, strict = TRUE))

  # TODO: What if it is zero?

  # Check that the maximum number of unobserved species is a positive number
  testit::assert(is_positive_integer(nmax))

  # Check the tolerance criterion
  testit::assert(is_positive(tol))

  # We must turn the parameters into a vector
  pars <- unlist(pars)

  # Transformation functions
  trans <- range_transform
  untrans <- function(x) range_transform(x, inverse = TRUE)

  # Re-scale the parameters
  pars <- trans(pars)

  # TODO: Make the transformation optional? Or at least let the user choose
  # wether they want the simplex to untransform.

  # Extra arguments (names as expected by the likelihood function)
  extra <- list(
    data = data, island_age = island_age, M = M, nmax = nmax,
    control = control_ode
  )

  # TODO: Add verbose to the optimizer, or choose where to have it at least.

  # Optimize
  optimizer(calc_loglik, pars, control_ml, method, extra)

}
