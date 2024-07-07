# Check parameters passed to a function
check_pars <- function(pars) {

  # pars: the parameters

  # Must be a list
  if (!is.list(pars)) stop("pars must be a list")

  # Parameters must be named
  if (is.null(names(pars))) stop("pars must be a named list")

  # Parameter names
  par_names <- c("lambda_c", "mu", "gamma", "lambda_a")

  # Find mismatching items
  is_unknown <- !(names(pars) %in% par_names)
  is_missing <- !(par_names %in% names(pars))

  # Error if any
  if (any(is_unknown)) stop("Unknown parameter(s): ", paste(names(pars)[is_unknown], collapse = ", "))
  if (any(is_missing)) stop("Missing parameter(s): ", paste(par_names[is_missing], collapse = ", "))

  # Check that all parameters are positive numbers
  if (!is_positive_vector(unlist(pars))) stop("pars must be a list of positive numbers")

}

# Function to check that the timings are good
check_times <- function(tcol, tmin, tmax, branching_times, island_age, clade) {

  # tcol: known colonization time (could be NULL)
  # tmin: minimum colonization time (could be NULL)
  # tmax: maximum colonization time (could be NULL)
  # branching_times: vector of branching times (could be empty)
  # island_age: age of the island
  # clade: clade name

  # Check island age
  testit::assert(is_negative(island_age))
  testit::assert(is.character(clade))

  # Error message preamble
  in_clade <- paste("in clade", clade)

  # There must be colonization otherwise there is no clade
  if (is.null(tcol) & is.null(tmax))
    stop(in_clade, ": tcol and tmax cannot both be NULL")

  # If known colonization time, upper and lower bounds should not be set
  if (!is.null(tcol) & !is.null(tmax))
    stop(in_clade, ": if tcol is provided then tmax must be NULL")

  # If a lower bound is given, then an upper bound should have been given too
  if (!is.null(tmin) & is.null(tmax))
    stop(in_clade, ": if tmin is provided then tmax must be too")

  # Check that the colonization times are all negative numbers
  if (!is.null(tcol) & !is_negative(tcol)) stop(in_clade, ": tcol must be a negative number")
  if (!is.null(tmax) & !is_negative(tmax)) stop(in_clade, ": tmax must be a negative number")
  if (!is.null(tmin) & !is_negative(tmin)) stop(in_clade, ": tmin must be a negative number")

  # Is there any branching time?
  is_branching <- length(branching_times) > 0L

  # And that if provided branching times are too
  if (is_branching & !all(is_negative_vector(branching_times)))
    stop(in_clade, ": branching_times must all be negative numbers")

  # Branching times should have been sorted already
  if (is_branching) testit::assert(all(diff(branching_times) > 0))

  # TODO: What happens when some branching times are simultaneous?

  # The first (oldest) branching time
  bt1 <- branching_times[1]

  # Check order of timings
  if (!is.null(tmax)) if (island_age > tmax) stop(in_clade, ": tmax must be younger than island_age")
  if (!is.null(tcol)) if (island_age > tcol) stop(in_clade, ": tcol must be younger than island_age")
  if (!is.null(tmin)) if (tmax > tmin) stop(in_clade, ": tmin must be younger than tmax")
  if (is_branching) if (island_age > bt1) stop(in_clade, ": branching_times must be younger than island_age")
  if (!is.null(tmax) & is_branching) if (tmax > bt1) stop(in_clade, ": branching_times must be younger than tmax")
  if (!is.null(tcol) & is_branching) if (tcol > bt1) stop(in_clade, ": branching_times must be younger than tcol")
  if (!is.null(tmin) & is_branching) if (tmin > bt1) stop(in_clade, ": branching_times must be younger than tmin")

  # TODO: Special case of is_number() for non-infinites?

  # TODO: What happens when some of those times are simultaneous?

}

# Function to check a clade
check_clade <- function(clade, name, island_age) {

  # clade: the clade
  # name: name of the clade
  # island_age: age of the island

  # Checks
  testit::assert(is_negative(island_age))
  testit::assert(is.character(name))

  # Error message preamble
  clade_name <- paste("clade", name)

  # Must be a list
  if (!is.list(clade)) stop(clade_name, " must be a list")

  # Must contain five elements
  if (length(clade) == 0L) stop(clade_name, " cannot be empty")

  # Items must have names
  if (is.null(names(clade))) stop(clade_name, " must be a named list")

  # Must contain the following items
  clade_items <- c("is_present", "tcol", "branching_times", "tmax", "tmin")

  # Find mismatching items
  is_unknown <- !(names(clade) %in% clade_items)
  is_missing <- !(clade_items %in% names(clade))

  # Error if any
  if (any(is_unknown)) stop(clade_name, " has unknown element(s): ", paste(names(clade)[is_unknown], collapse = ", "))
  if (any(is_missing)) stop(clade_name, " has missing element(s): ", paste(clade_items[is_missing], collapse = ", "))

  # Make sure whether the mainland colonist is present is a yes or no
  if (!is_yes_no(clade$is_present)) stop("in ", clade_name, ": is_present must be TRUE or FALSE")

  # Check the timings
  with(clade, check_times(tcol, tmin, tmax, branching_times, island_age, name))

}

# Function to check the input data
check_data <- function(data, island_age) {

  # data: the list of clades
  # island_age: the island age

  # Check
  testit::assert(is_negative(island_age))

  # Must be a list
  if (!is.list(data)) stop("data must be a list")

  # Must have at least one colonist clade
  if (length(data) == 0L) stop("data must at least have one element")

  # TODO: What about an empty island?

  # Clades must be named
  if (is.null(names(data))) stop("data must be a named list")

  # Each element in the list must be a clade
  for (i in seq(data)) check_clade(clade = data[[i]], name = names(data)[[i]], island_age)

}

# Function to find maximum likelihood estimates
daisie_ml <- function(

  data, pars, island_age, M, nmax,
  control_ml = list(),
  control_ode = list(),
  method = "subplex",
  verbose = FALSE

) {

  # data: the input data
  # pars: model parameters
  # island_age: age of the island
  # M: size of the mainland pool
  # nmax: maximum allowed number of unobserved species
  # control_ml: named list of options for the likelihood maximization algorithm
  # control_ode: named list of options for the integrator
  # method: which method ("simplex" or "subplex") to use for optimization?
  # verbose: whether to display messages

  # TODO: Implement correction for clades older than island age and sorting
  # of branching times, but those might be done in the preparation function.

  # Checks
  if (!is_negative(island_age)) stop("island_age must be negative")
  if (!is_positive_integer(M, strict = TRUE)) stop("M must be a positive integer")
  if (!is_positive_integer(nmax)) stop("nmax must be a positive integer")
  if (!is.character(method)) stop("method must be a character string")
  if (!is_yes_no(verbose)) stop("verbose must be TRUE or FALSE")

  # Check some likelihood maximization parameters
  if (!is.list(control_ml)) stop("control_ml must be a list")
  if ("trans" %in% names(control_ml)) stop("trans cannot be provided in control_ml")
  if ("untrans" %in% names(control_ml)) stop("untrans cannot be provided in control_ml")

  # Note: the rest will be checked by the optimizer.

  # TODO: Figure where to turn island age into a positive number.

  # TODO: What if M is zero?

  # TODO: Where do we check for ODE solver control options?

  # TODO: Get our shit together with where we check for what.

  # TODO: Think about data preparation.

  # TODO: Add diversity-dependence.

  # Check the input data
  check_data(data, island_age)

  # Check the parameters
  check_pars(pars)

  # We must turn the parameters into a vector
  pars <- unlist(pars)

  # Transformation functions
  trans <- range_transform
  untrans <- function(x) range_transform(x, inverse = TRUE)

  # Add those functions to the control options
  control_ml <- c(control_ml, trans = trans, untrans = untrans)

  # Re-scale the parameters
  pars <- trans(pars)

  # TODO: Make the transformation optional? Or at least let the user choose
  # wether they want the simplex to untransform.

  # Extra arguments (names as expected by the likelihood function)
  extra <- list(
    data = data, island_age = island_age, M = M, nmax = nmax,
    control = control_ode
  )

  # Optimize
  optimizer(calc_loglik, pars, control_ml, extra, method, verbose, warn = FALSE)

}

# TODO: Could we not just use subplex?
