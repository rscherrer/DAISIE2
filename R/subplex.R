# Function to set the control options
make_control_subplex <- function(options) {

  # options: named list of user-specified options

  # Initial checks
  testit::assert(is.list(options))

  # Default control parameters
  control <- list(

    # Numerical tolerance parameters to check convergence
    abstol = 1e-6,
    reltol = 1e-4,

    # Maximum number of iterations
    maxit = 1000L,

    # How much to dodge parameters equal to one half?
    jitter = 0

  )

  # Return default parameters if nothing specified
  if (length(options) == 0L) return(control)

  # Check names in the input
  testit::assert(!is.null(names(options)))
  testit::assert(names(options) %in% names(control))

  # Update non-default parameters
  control[names(options)] <- options

  return(control)

}

# Function to check the control parameters
check_control_subplex <- function(control) {

  # control: list containing the control parameters

  # Check that it is a list
  testit::assert(is.list(control))

  # Names of the control parameters
  control_names <- c("abstol", "reltol", "maxit", "jitter")

  # Check that they are all here
  testit::assert(all(names(control) %in% control_names))
  testit::assert(all(control_names %in% names(control)))

  # With the control parameters...
  with(control, {

    # Perform checks
    testit::assert(is_positive(abstol))
    testit::assert(is_positive(reltol))
    testit::assert(is_positive_integer(maxit, strict = TRUE))
    testit::assert(is_number(jitter))

  })
}

# Wrapper around a function that implements subplex optimization
subplex <- function(fun, pars, control = list(), extra = list()) {

  # fun: the function to maximize
  # pars: vector of parameters to optimize (initial guesses, can be named)
  # control: a list of named options for the algorithm, as defined in control()
  # extra: named list of extra arguments of the function to optimize

  # Update default options with user choices
  control <- make_control_subplex(control)

  # Check the control parameters
  check_control_subplex(control)

  # Dodge some values if needed
  pars[pars == 0.5] <- 0.5 - control$jitter

  # Prepare the function to minimize (in form expected by subplex)
  this_fun <- function(pars, extra) { -call_fun(fun, pars, extra) }

  # Pass all that to the subplex routine
  out <- subplex::subplex(

    par = pars,
    fn = this_fun,
    control = control,
    extra = extra

  )

  # TODO: Warnings were suppressed in original code.

  # Combine the relevant output in a list
  return(with(out, list(pars = par, fvalue = -value, conv = convergence)))

}
