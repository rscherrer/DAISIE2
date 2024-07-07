# Perform argument checks
check_args <- function(fun, pars, control, extra, verbose) {

  # Check
  testit::assert(is.function(fun))
  testit::assert(is_numeric_vector(pars))
  testit::assert(is.list(control))
  testit::assert(length(control) == 0L | !is.null(names(control)))
  testit::assert(is.list(extra))
  testit::assert(is_yes_no(verbose))

}

# Default options for simplex
get_default_options_simplex <- function() {

  list(

    # Numerical tolerance parameters to check convergence
    rtol = 1e-4,
    ftol = 1e-5,
    atol = 1e-6,

    # Maximum number of iterations
    maxiter = 1000L,

    # By how much to shift the initial guesses to initialize the simplex (e.g. 0.05 for 5%)
    delta = 0.05,

    # Value by which to shift in case initial guess is zero
    dzero = 0.00025,

    # Geometric transformation parameters
    rho = 1,
    chi = 2,
    psi = 0.5,
    sigma = 0.5,

    # Transformation functions used to (back-)transform the data
    trans = NULL,
    untrans = NULL

  )
}

# Default options for subplex
get_default_options_subplex <- function() {

  list(

    # Numerical tolerance parameters to check convergence
    atol = 1e-6,
    rtol = 1e-4,

    # Maximum number of iterations
    maxiter = 1000L,

    # How much to dodge parameters equal to one half?
    jitter = 0,

    # Function to untransform (re-scaled) parameters if needed
    untrans = NULL

  )
}

# Default options for the optimizer
get_default_options_optimizer <- function() {

  list(

    # Number of optimization cycles
    ncycles = 1L,

    # Numerical tolerance to check convergence between cycles
    ctol = 1e-6

  )
}

# Function to set the control options
make_control <- function(options, method, meta = FALSE, verbose = FALSE) {

  # options: named list of user-specified options
  # method: the algorithm to use
  # meta: whether to append general optimizer options
  # verbose: whether to display messages

  # Note: calling the algorithm from within the optimizer allows to run it
  # for multiple cycles, but the optimizer needs its own extra options.

  # Initial check
  testit::assert(is.list(options))
  testit::assert(is.character(method))
  testit::assert(method %in% c("simplex", "subplex"))
  testit::assert(is_yes_no(meta))
  testit::assert(is_yes_no(verbose))

  # Default control parameters
  control <- switch(
    method,
    simplex = get_default_options_simplex(),
    subplex = get_default_options_subplex()
  )

  # Add meta-options if needed
  if (meta) control <- c(control, get_default_options_optimizer())

  # Return default parameters if nothing specified
  if (length(options) == 0L) return(control)

  # Check names in the input
  testit::assert(!is.null(names(options)))

  # Are there any unknown options?
  is_unknown <- !(names(options) %in% names(control))

  # If needed...
  if (verbose & any(is_unknown)) {

    # Identify them
    unknown_options <- names(options)[is_unknown]

    # Tell the user (those may be typos)
    message("Unknown option(s) (will be ignored): ", paste(unknown_options, collapse = " "))

  }

  # If no known options are specified, exit
  if (sum(is_unknown) == length(is_unknown)) return(control)

  # Update non-default parameters
  control[names(options)] <- options

  return(control)

}

# Function to check the transformation functions
check_trans <- function(trans, untrans) {

  # trans: the transformation function
  # untrans: its inverse

  # Transformation functions must be either functions or nothing
  if (!is.function(trans) & !is.null(trans)) stop("trans must be NULL or a function")
  if (!is.function(untrans) & !is.null(untrans)) stop("untrans must be NULL or a function")

  # If one transformation is provided then both must be
  if (is.function(trans) != is.function(untrans))
    stop("both trans and untrans must be provided, or none")

  # Exit now if not provided
  if (is.null(trans)) return()

  # Dummy value
  dummy <- 0.9

  # They must be inverse of each other (try with a dummy value)
  if (trans(untrans(dummy)) != dummy) stop("untrans is not the inverse function of trans")

}

# Function to check the control parameters
check_control_simplex <- function(control) {

  # control: list containing the control parameters

  # Check that it is a list
  testit::assert(is.list(control))

  # Names of the control parameters
  control_names <- c(
    "rtol", "ftol", "atol", "maxiter", "delta", "dzero",
    "rho", "chi", "psi", "sigma", "trans", "untrans"
  )

  # Check that they are all here
  testit::assert(all(control_names %in% names(control)))

  # With the control parameters...
  with(control, {

    # Perform checks
    if (!is_positive(rtol)) stop("rtol must be a positive number")
    if (!is_positive(ftol)) stop("ftol must be a positive number")
    if (!is_positive(atol)) stop("atol must be a positive number")
    if (!is_positive_integer(maxiter)) stop("maxiter must be a positive integer")
    if (!is_positive(delta)) stop("delta must be a positive number")
    if (!is_number(dzero)) stop("dzero must be a number")
    if (!is_positive(rho)) stop("rho must be a positive number")
    if (!is_positive(chi)) stop("chi must be a positive number")
    if (!is_positive(psi)) stop("psi must be a positive number")
    if (!is_positive(sigma)) stop("sigma must be a positive number")

    # Check the transformation functions
    check_trans(trans, untrans)

  })
}

# Function to assemble the output of the optimizer in a list
prepare_output <- function(pars, fvalue, conv, untrans) {

  # pars: the vector of parameters
  # fvalue: function value for those parameters
  # conv: convergence code (0 if yes, -1 if iterations exceeded, plus others)
  # untrans: function to untransform the parameters if re-scaled

  # Note: convergence can be zero or minus one for both simplex and subplex
  # but subplex offers other non-convergence codes (see ?subplex::subplex).

  list(
    pars = if (is.null(untrans)) pars else untrans(pars),
    fvalue = fvalue,
    conv = conv
  )
}

# Function to return a list of information for a new vertex
new_vertex <- function(vertex, fvalue, how) {

  # vertex: the coordinates of the vertex
  # fvalue: the function value at that vertex
  # how: what geometric transformation does that vertex correspond to?

  # Combine the info into a list
  return(list(vertex = vertex, fvalue = fvalue, how = how))

}

# Check a function value and error if needed
check_fvalue <- function(f) {

  # Error if needed
  if (is.na(f)) stop("Problem encountered in function evaluation.")

  # TODO: Any other actions in case of NA? Like early exit? With -Inf log-lik.?

}

# Function to update the worst vertex in the simplex
get_updated_vertex <- function(

  x, xbar, fun, fmin, fstm, fmax, rho, chi, psi, sigma,
  untrans = NULL, extra = list()

) {

  # x: coordinates of the worst vertex
  # xbar: centroid of the edge to flip along
  # fun: the function to measure
  # fmin: minimum (i.e. best) score out of all vertices
  # fstm: second-to-maximum score out of all vertices
  # fmax: maximum (i.e. worst) score out of all vertices
  # rho, chi, psi, sigma: update parameters
  # untrans: the function to untransform (re-scaled) parameters if needed
  # extra: named list of extra arguments for the function to be optimized

  # Checks that were not already performed
  testit::assert(is_number(x, scalar = FALSE))
  testit::assert(is_number(xbar, scalar = FALSE))
  testit::assert(length(x) == length(xbar))

  # Vector of reference function values
  fvalues <- c(fmin, fstm, fmax)

  # Check that those values are in increasing order
  testit::assert(all(is_number(fvalues, scalar = FALSE)))
  testit::assert(all(diff(fvalues) >= 0))

  # Find the reflection of the worst point wrt the centroid
  xr <- (1 + rho) * xbar - rho * x
  fxr <- -call_fun(fun, xr, extra, untrans)
  check_fvalue(fxr)

  # TODO: Does that mean we can make it work with lists of parameters?

  # Set reflection as our new prospect
  out <- new_vertex(xr, fxr, how = "reflect")

  # If the reflection is a new minimum...
  if (fxr < fmin) {

    # Expand it to see if we can minimize further
    xe <- (1 + rho * chi) * xbar - rho * chi * x
    fxe <- -call_fun(fun, xe, extra, untrans)
    check_fvalue(fxe)

    # Return the expansion if smaller
    if (fxe < fxr) return(new_vertex(xe, fxe, how = "expand"))

    # Otherwise stick with the reflection
    return(out)

  }

  # Else, stick with the reflection if it is better than the worst two
  if (fxr < fstm) return(out)

  # Else, if it is at least not the worst...
  if (fxr < fmax) {

    # Contract the worst vertex
    xco <- (1 + psi * rho) * xbar - psi * rho * x
    fxco <- -call_fun(fun, xco, extra, untrans)
    check_fvalue(fxco)

    # Keep the contraction if it does better
    if (fxco <= fxr) return(new_vertex(xco, fxco, how = "contract_outside"))

    # Otherwise shrink the simplex
    return(new_vertex(NULL, NULL, how = "shrink"))

  }

  # If the reflection was the worst, then compute the inside contraction
  xci <- (1 - psi) * xbar + psi * x
  fxci <- -call_fun(fun, xci, extra, untrans)
  check_fvalue(fxci)

  # Keep this contraction if it does better
  if (fxci < fmax) return(new_vertex(xci, fxci, how = "contract_inside"))

  # Otherwise shrink the simplex
  return(new_vertex(NULL, NULL, "shrink"))

}

# Function to compute (Nelder-Mead) simplex optimization
simplex <- function(fun, pars, control = list(), extra = list(), verbose = FALSE) {

  # fun: the function to maximize
  # pars: vector of parameters to optimize (initial guesses, can be named)
  # control: a list of named options for the algorithm, as defined in control()
  # extra: (named) list of extra arguments of the function to optimize
  # verbose: whether to display messages

  # Note: here we assume that the function to optimize takes a vector of
  # parameters as one of its arguments ("pars"), and those are the parameters
  # to optimize. The extra parameters in "extra" should be named so they can
  # be matched to the right argument in the function. Once those are passed,
  # the parameters to optimize will be passed to the first argument not yet
  # passed in the function to optimize. If the arguments in "extra" are not
  # named, they will be passed in order to the function.

  # Check arguments in common with other optimizer functions
  check_args(fun, pars, control, extra, verbose)

  # Update default options with user choices
  control <- make_control(control, method = "simplex", verbose = verbose)

  # Check the control parameters
  check_control_simplex(control)

  # Number of parameters (i.e. dimensions)
  npars <- length(pars)

  # Number of vertices
  nvertices <- npars + 1L

  # Prepare the simplex as a matrix of coordinates for each vertex
  V <- matrix(rep(pars, nvertices), nrow = npars)

  # Ratio by which to dodge each vertex dimension
  r0 <- 1 + control$delta

  # Make a copy of the parameters
  dpars <- pars

  # Extract the transformation functions
  trans <- control$trans
  untrans <- control$untrans

  # Untransform if needed
  if (!is.null(untrans)) dpars <- untrans(pars)

  # Dodge
  dpars <- dpars * r0

  # Transform them back if needed
  if (!is.null(trans)) dpars <- trans(dpars)

  # Prepare a vector of dodge factors
  r <- rep(1, npars)

  # That ratio cannot be measured on zero-valued parameters so identify those
  inonzero <- pars != 0

  # Measure the dodge ratios on the transformed scale
  r[inonzero] <- dpars[inonzero] / pars[inonzero]

  # Should have run smoothly
  testit::assert(!any(is.nan(r)))

  # Make sure that they are no smaller than the ratio on the untransformed scale
  r[r < r0] <- r0

  # Indices of the coordinates to add a little deviation to
  ii <- seq(npars)
  ii <- ii * npars + ii

  # Note: those indices are positions in the matrix corresponding to (i, i + 1)
  # for as many i's as there are dimensions (parameters). This way, each
  # vertex is shifted by a small amount in only one particular direction.
  # See https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method for
  # a description of the algorithm.

  # Shift each vertex by the specified ratio along only one dimension
  V[ii] <- r * V[ii]

  # Small deviation for zero-valued coordinates
  V[ii][V[ii] == 0] <- control$dzero

  # Add parameter names to the simplex if needed
  if (!is.null(names(pars))) rownames(V) <- names(pars)

  # Compute the function for each vertex
  fvalues <- apply(V, 2L, \(x) { -call_fun(fun, x, extra, untrans) })

  # Note: this algorithm maximizes a function by finding the minimum
  # of its negative.

  # Check that all could be computed
  for (f in fvalues) check_fvalue(f)

  # Reorder vertices from minimum to maximum function value
  jj <- order(fvalues)
  V <- V[, jj]
  if (npars == 1L) V <- matrix(V, nrow = npars)
  fvalues <- fvalues[jj]

  # Note: this way the last vertex is the worst one.

  # Initialize counter
  iter <- 1L

  # For each iteration...
  while (iter <= control$maxiter) {

    # Compute useful metrics for convergence checking
    maxdf <- max(abs(fvalues - fvalues[1]))
    dV <- abs(V - V[, 1])

    # Note: this latter line compute the distance between each vertex (column)
    # and the first one.

    # Measure convergence in both function evaluation and simplex coordinates
    is_conf <- !is.nan(maxdf) & maxdf <= control$ftol * abs(fvalues[1])
    is_conv <- with(control, max(dV - rtol * abs(V[, 1])) <= 0 & max(dV) <= atol)

    # Break the loop if convergence is reached
    if (is_conf | is_conv) break

    # Centroid of the n (not n + 1)-dimensional edge along which to flip the simplex
    xbar <- rowSums(as.matrix(V[, 1:npars])) / npars

    # Note: this would be a line in a 2-dimensional parameter space here
    # the simplex is a triangle, for example, but it would be a triangle
    # in a 3-dimensional parameter space where the simplex is a tetrahedron.

    # Note: we turn the object into a matrix again because rowSums errors
    # when only one parameter (and it is applied on a scalar).

    # Figure how to update the worst vertex
    new_vertex <- with(control, get_updated_vertex(

      # Vertex to update and centroid wrt which to transform
      V[, nvertices], xbar,

      # The function
      fun,

      # Function values to compare to
      fmin = fvalues[1], fstm = fvalues[npars], fmax = fvalues[nvertices],

      # Geometric transformation parameters
      rho, chi, psi, sigma,

      # Function to untransform (re-scaled) parameters if needed
      untrans,

      # Extra arguments for the function to optimize
      extra

    ))

    # If the simplex must be shrunk...
    if (new_vertex$how == "shrink") {

      # For each vertex except the best one
      kk <- 2:nvertices

      # Shrink the simplex by shifting that vertex towards the best one
      V[, kk] <- V[, 1] + control$sigma * (V[, kk] - V[, 1])
      fvalues[kk] = apply(as.matrix(V[, kk]), 2L, \(x) -call_fun(fun, x, extra, untrans))

      # Note: we force into a matrix in case there is only one element.

      # Check function values for problems
      for (f in fvalues) check_fvalue(f)

      # Check parameter values for problems
      testit::assert(all(apply(as.matrix(V[, kk]), 2L, is_numeric_vector)))

    } else {

      # Otherwise update the worst vertex
      V[, nvertices] <- new_vertex$vertex
      fvalues[nvertices] <- new_vertex$fvalue

      # Check parameter values for problems
      testit::assert(is_numeric_vector(V[, nvertices]))

    }

    # Reorder vertices from minimum to maximum function value
    jj <- order(fvalues)
    V <- V[, jj]
    if (npars == 1L) V <- matrix(V, nrow = npars)
    fvalues <- fvalues[jj]

    # Verbose if needed
    if (verbose) {

      # Extract elements to display
      best_pars <- if (is.null(untrans)) V[, 1] else untrans(V[, 1])
      best_fvalue <- -fvalues[1]

      # Display
      message(
        "Iter. ", iter,
        ", oper.: ", new_vertex$how,
        ", param.: ", paste(best_pars, collapse = ", "),
        ", value: ", best_fvalue
      )
    }

    # Stop if we reached too large negative numbers
    if (any(fvalues == -Inf)) break

    # Increment counter
    iter <- iter + 1L

  }

  # Verbose if needed
  if (verbose & iter <= control$maxiter) message("Optimization terminated successfully.")
  if (verbose & iter > control$maxiter) message("Maximum number of iterations exceeded")

  # Return output
  return(prepare_output(V[, 1], -fvalues[1], -as.integer(iter > control$maxiter), untrans))

}

# Function to check the control parameters
check_control_subplex <- function(control) {

  # control: list containing the control parameters

  # Check that it is a list
  testit::assert(is.list(control))

  # Names of the control parameters
  control_names <- c("atol", "rtol", "maxiter", "jitter", "untrans")

  # Check that they are all here
  testit::assert(all(control_names %in% names(control)))

  # With the control parameters...
  with(control, {

    # Perform checks
    if (!is_positive(atol)) stop("atol must be a positive number")
    if (!is_positive(rtol)) stop("rtol must be a positive number")
    if (!is_positive_integer(maxiter)) stop("maxiter must be a positive integer")
    if (!is_number(jitter)) stop("jitter must be a number")
    if (!is.null(untrans) & !is.function(untrans)) stop("untrans must be NULL or a function")

  })
}

# Wrapper around a function that implements subplex optimization
subplex <- function(fun, pars, control = list(), extra = list(), verbose = FALSE) {

  # fun: the function to maximize
  # pars: vector of parameters to optimize (initial guesses, can be named)
  # control: a list of named options for the algorithm, as defined in control()
  # extra: named list of extra arguments of the function to optimize
  # verbose: whether to display messages

  # Check arguments in common with other optimizer functions
  check_args(fun, pars, control, extra, verbose)

  # Update default options with user choices
  control <- make_control(control, method = "subplex", verbose = verbose)

  # Check the control parameters
  check_control_subplex(control)

  # Dodge some values if needed
  pars[pars == 0.5] <- 0.5 - control$jitter

  # Unpack the transformation function for re-scaled parameters if needed
  untrans <- control$untrans

  # Prepare the function to minimize (in form expected by subplex)
  this_fun <- function(pars, extra, untrans) { -call_fun(fun, pars, extra, untrans) }

  # If the optimization is not even to be started...
  if (control$maxiter == 0L) {

    # Compute function value
    fvalue <- this_fun(pars, extra, control$untrans)
    check_fvalue(fvalue)

    # Early exit with convergence code for max. iterations elapsed
    return(prepare_output(pars, -fvalue, conv = -1, untrans = control$untrans))

  }

  # Pass all that to the subplex routine
  out <- subplex::subplex(

    par = pars,
    fn = this_fun,
    control = with(control, list(reltol = rtol, abstol = atol, maxit = maxiter)),
    extra = extra,
    untrans = untrans

    # Note: the control options have to be renamed for subplex::subplex

  )

  # Combine the relevant output in a list
  return(with(out, prepare_output(par, -value, convergence, untrans)))

}

# Function to check the control parameters
check_control_optimizer <- function(control) {

  # control: list containing the control parameters

  # Check that it is a list
  testit::assert(is.list(control))

  # Names of the control parameters
  control_names <- c("ncycles", "ctol", "untrans")

  # Check that they are all here
  testit::assert(all(control_names %in% names(control)))

  # With the provided options...
  with(control, {

    # Check that they are fine
    if (!is_positive_integer(ncycles)) stop("ncycles must be a positive integer")
    if (!is_positive(ctol)) stop("ctol must be a positive number")
    if (!is.null(untrans) & !is.function(untrans)) stop("untrans must be NULL or a function")

  })

}

# Routine to find the maximum of a function
optimizer <- function(

  fun, pars, control = list(), extra = list(),
  method = "subplex", verbose = FALSE, warn = TRUE

) {

  # fun: the function to maximize
  # pars: vector of parameters to optimize (initial guesses, can be named)
  # control: a list of named options for the algorithm, as defined in control()
  # extra: named list of extra arguments of the function to optimize
  # method: which method (of "simplex" or "subplex") to use for optimization?
  # verbose: whether to display messages
  # warn: whether to display warnings from the underlying algorithm

  # Checks
  if (!is.function(fun)) stop("fun must be a function")
  if (!is_numeric_vector(pars)) stop("pars must be a numeric vector")
  if (!is.list(control)) stop("control must be a named list")
  if (length(control) > 0L & is.null(names(control))) stop("control must be a named list")
  if (!is.list(extra)) stop("extra must be a list")
  if (!is_yes_no(verbose)) stop("verbose must be TRUE or FALSE")

  # Methods allowed
  methods <- c("simplex", "subplex")

  # Check method
  if (!is.character(method)) stop("method must be ", paste(methods, collapse = " or "))
  if (!(method %in% methods)) stop("method must be ", paste(methods, collapse = " or "))

  # Update default options with user choices
  control <- make_control(control, method, meta = TRUE, verbose = verbose)

  # Check the control parameters
  check_control_optimizer(control)

  # Remove meta-parameters from the options being passed to the algorithm
  control0 <- control[!(names(control) %in% c("ncycles", "ctol"))]

  # Note: this is so the inner function does not complain about unknown
  # parameters.

  # Set the optimizer
  this_optimizer <- this_optimizer0 <- get(method)

  # Note: the transformation functions for re-scaled parameters are in
  # the control options because we want the different algorithms (subplex,
  # simplex) to take the same number of arguments so they can be called using
  # get() easily (and e.g. simplex takes one more transformation function
  # than subplex).

  # Turn off warnings if needed
  if (!warn) this_optimizer <- function(...) suppressWarnings(this_optimizer0(...))

  # Initialize counter
  cycle <- 1L

  # Initialize function value (to check convergence)
  lastf <- Inf

  # Compute function value in case of zero cycle
  fvalue <- call_fun(fun, pars, extra, control$untrans)
  check_fvalue(fvalue)

  # Output to return in case of zero cycle
  out <- prepare_output(pars, -fvalue, conv = -1, control$untrans)

  # For each optimization cycle...
  while (cycle <= control$ncycles) {

    # Verbose if needed
    if (verbose) message("Cycle ", cycle)

    # Optimize the likelihood function
    out <- this_optimizer(fun = fun, pars = pars, control = control0, extra = extra, verbose = verbose)

    # Use optimized parameters as starting point for the next cycle
    pars <- out$pars
    newf <- out$fvalue

    # Check function value for error
    check_fvalue(newf)

    # Exit the loop if function values are not comparable
    if (is.nan(newf - lastf)) break

    # Note: e.g. if they are both infinite and of the same sign.

    # Compute convergence criterion between cycles
    has_conv <- abs(newf - lastf) <= control$ctol

    # Break the loop if converged
    if (has_conv) break

    # Update the latest function value
    lastf <- newf

    # Increment cycles
    cycle <- cycle + 1L

  }

  # Non-convergence criterion
  non_conv <- verbose & cycle > control$ncycles & out$conv != 0L

  # Verbose if needed
  if (verbose)
    message(ifelse(non_conv, "More cycles in optimization recommended.", "No more cycles needed."))

  return(out)

}

# TODO: Verify punctuation in messages at the end.
