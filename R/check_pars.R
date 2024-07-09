# Check parameter values
check_par_values <- function(pars) {

  # Check that...
  with(pars, {

    # ... all parameters are positive numbers
    if (!is_positive(lambda_c)) stop("lambda_c must be positive")
    if (!is_positive(mu)) stop("mu must be positive")
    if (!is_positive(gamma)) stop("gamma must be positive")
    if (!is_positive(lambda_a)) stop("lambda_a must be positive")

  })

  # TODO: What if they are zero?

}

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

  # Check parameter values
  check_par_values(pars)

}
