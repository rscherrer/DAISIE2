# Default options for simplex
get_default_options_simplex <- function() {

  list(

    # Numerical tolerance parameters to check convergence
    rtol = 1e-4,
    ftol = 1e-5,
    atol = 1e-6,

    # Maximum number of iterations
    maxiter = 1000L,

    # TODO: Previously 1000 * round((1.25)^npars for simplex

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
make_control <- function(options, method, meta = FALSE) {

  # options: named list of user-specified options
  # method: the algorithm to use
  # meta: whether to append general optimizer options

  # Note: calling the algorithm from within the optimizer allows to run it
  # for multiple cycles, but the optimizer needs its own extra options.

  # Initial check
  testit::assert(is.list(options))
  testit::assert(is.character(method))
  testit::assert(method %in% c("simplex", "subplex"))
  testit::assert(is_yes_no(meta))

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

  # Or if no known options are specified
  if (!any(names(options) %in% names(control))) return(control)

  # TODO: Maybe issue a message or warning when verbose to say that some
  # control option was ignored (same in simplex and subplex).

  # Update non-default parameters
  control[names(options)] <- options

  return(control)

}

# Function to check the transformation functions
check_trans <- function(trans, untrans) {

  # trans: the transformation function
  # untrans: its inverse

  # Transformation functions must be either functions or nothing
  testit::assert(is.function(trans) | is.null(trans))
  testit::assert(is.function(untrans) | is.null(untrans))

  # If one transformation is provided then both must be
  testit::assert(is.function(trans) == is.function(untrans))

  # Exit now if not provided
  if (is.null(trans)) return()

  # Dummy value
  dummy <- 0.9

  # They must be inverse of each other (try with a dummy value)
  testit::assert(trans(untrans(dummy)) == dummy)

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
    testit::assert(is_number(rtol, sign = 1))
    testit::assert(is_number(ftol, sign = 1))
    testit::assert(is_number(atol, sign = 1))
    testit::assert(is_number(maxiter, integer = TRUE, sign = 1))
    testit::assert(is_number(delta, sign = 1))
    testit::assert(is_number(dzero))
    testit::assert(all(is_number(c(rho, chi, psi, sigma), scalar = FALSE, sign = 1)))

    # Check the transformation functions
    check_trans(trans, untrans)

  })
}

# Function to return a list of information for a new vertex
new_vertex <- function(vertex, fvalue, how) {

  # vertex: the coordinates of the vertex
  # fvalue: the function value at that vertex
  # how: what geometric transformation does that vertex correspond to?

  # Combine the info into a list
  return(list(vertex = vertex, fvalue = fvalue, how = how))

}

# Function to update the worst vertex in the simplex
get_updated_vertex <- function(x, xbar, fun, fmin, fstm, fmax, rho, chi, psi, sigma, untrans = NULL, extra = list()) {

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

  # TODO: Does that mean we can make it work with lists of parameters?

  # Set reflection as our new prospect
  out <- new_vertex(xr, fxr, how = "reflect")

  # TODO: Sometimes the function values are not numbers, what do we do then?

  # If the reflection is a new minimum...
  if (fxr < fmin) {

    # Expand it to see if we can minimize further
    xe <- (1 + rho * chi) * xbar - rho * chi * x
    fxe <- -call_fun(fun, xe, extra, untrans)

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

    # Keep the contraction if it does better
    if (fxco <= fxr) return(new_vertex(xco, fxco, how = "contract_outside"))

    # Otherwise shrink the simplex
    return(new_vertex(NULL, NULL, how = "shrink"))

  }

  # If the reflection was the worst, then compute the inside contraction
  xci <- (1 - psi) * xbar + psi * x
  fxci <- -call_fun(fun, xci, extra, untrans)

  # Keep this contraction if it does better
  if (fxci < fmax) return(new_vertex(xci, fxci, how = "contract_inside"))

  # Otherwise shrink the simplex
  return(new_vertex(NULL, NULL, "shrink"))

}

# Function to compute (Nelder-Mead) simplex optimization
simplex <- function(fun, pars, control = list(), extra = list()) {

  # fun: the function to maximize
  # pars: vector of parameters to optimize (initial guesses, can be named)
  # control: a list of named options for the algorithm, as defined in control()
  # extra: named list of extra arguments of the function to optimize

  # Note: here we assume that the function to optimize takes a vector of
  # parameters as one of its arguments ("pars"), and those are the parameters
  # to optimize. The extra parameters in "extra" should be named so they can
  # be matched to the right argument in the function. Once those are passed,
  # the parameters to optimize will be passed to the first argument not yet
  # passed in the function to optimize.

  # TODO: Specify starting values, stop if not.

  # TODO: Verbose.

  # Checks
  testit::assert(is.function(fun))
  testit::assert(is_number(pars, scalar = FALSE))

  # Number of parameters (i.e. dimensions)
  npars <- length(pars)

  # Update default options with user choices
  control <- make_control(control, method = "simplex")

  # Check the control parameters
  check_control_simplex(control)

  # Check that the extra arguments are a list, named if needed
  testit::assert(is.list(extra))
  if (length(extra) > 0L) testit::assert(!is.null(names(extra)))

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

  # TODO: [URGENT] We need to pass untransformed parameters to the function.

  # Compute the function for each vertex
  fvalues <- apply(V, 2L, \(x) { -call_fun(fun, x, extra, untrans) })

  # Note: this algorithm maximizes a function by finding the minimum
  # of its negative.

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

    } else {

      # Otherwise update the worst vertex
      V[, nvertices] <- new_vertex$vertex
      fvalues[nvertices] <- new_vertex$fvalue

    }

    # Reorder vertices from minimum to maximum function value
    jj <- order(fvalues)
    V <- V[, jj]
    if (npars == 1L) V <- matrix(V, nrow = npars)
    fvalues <- fvalues[jj]

    # Stop if we reached too large negative numbers
    if (any(fvalues == -Inf)) break

    # Increment counter
    iter <- iter + 1L

  }

  # Return...
  return(list(

    # The best parameters found and their associated function value
    pars = V[, 1],
    fvalue = -fvalues[1],

    # Whether convergence has been reached (code zero if yes)
    conv = as.integer(iter > control$maxiter)

  ))

  # TODO: Figure why many outputs are wrapped in invisible().

  # TODO: Make simplex() its own package? Or have a precise help file and export it.

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
    testit::assert(is_positive(atol))
    testit::assert(is_positive(rtol))
    testit::assert(is_positive_integer(maxiter, strict = TRUE))
    testit::assert(is_number(jitter))
    testit::assert(is.null(untrans) | is.function(untrans))

  })
}

# Wrapper around a function that implements subplex optimization
subplex <- function(fun, pars, control = list(), extra = list()) {

  # fun: the function to maximize
  # pars: vector of parameters to optimize (initial guesses, can be named)
  # control: a list of named options for the algorithm, as defined in control()
  # extra: named list of extra arguments of the function to optimize

  # Checks
  testit::assert(is.function(fun))
  testit::assert(is_number(pars, scalar = FALSE))

  # Check that the extra arguments are a list, named if needed
  testit::assert(is.list(extra))
  if (length(extra) > 0L) testit::assert(!is.null(names(extra)))

  # Update default options with user choices
  control <- make_control(control, method = "subplex")

  # Check the control parameters
  check_control_subplex(control)

  # Dodge some values if needed
  pars[pars == 0.5] <- 0.5 - control$jitter

  # Unpack the transformation function for re-scaled parameters if needed
  untrans <- control$untrans

  # Prepare the function to minimize (in form expected by subplex)
  this_fun <- function(pars, extra, untrans) { -call_fun(fun, pars, extra, untrans) }

  # Pass all that to the subplex routine
  out <- subplex::subplex(

    par = pars,
    fn = this_fun,
    control = with(control, list(reltol = rtol, abstol = atol, maxit = maxiter)),
    extra = extra,
    untrans = untrans

    # Note: the control options have to be renamed for subplex::subplex

  )

  # TODO: Warnings were suppressed in original code.

  # Combine the relevant output in a list
  return(with(out, list(pars = par, fvalue = -value, conv = convergence)))

}

# Function to check the control parameters
check_control_optimizer <- function(control) {

  # control: list containing the control parameters

  # Check that it is a list
  testit::assert(is.list(control))

  # Names of the control parameters
  control_names <- c("ncycles", "ctol")

  # Check that they are all here
  testit::assert(all(control_names %in% names(control)))

  # With the provided options...
  with(control, {

    # Check that they are fine
    testit::assert(is_positive_integer(ncycles, strict = TRUE))
    testit::assert(is_positive(ctol))

  })

  # TODO: What if a tolerance parameter is zero?

}

# Routine to find the maximum of a function
optimizer <- function(fun, pars, control = list(), method = "subplex", extra = list()) {

  # TODO: Use ellipsis for transformation functions?

  # fun: the function to maximize
  # pars: vector of parameters to optimize (initial guesses, can be named)
  # control: a list of named options for the algorithm, as defined in control()
  # method: which method (of "simplex" or "subplex") to use for optimization?
  # extra: named list of extra arguments of the function to optimize

  # Check the method supplied
  testit::assert(is.character(method))
  testit::assert(method %in% c("simplex", "subplex"))

  # Update default options with user choices
  control <- make_control(control, method, meta = TRUE)

  # Check the control parameters
  check_control_optimizer(control)

  # Set the optimizer
  this_optimizer <- get(method)

  # Initialize counter
  cycle <- 1L

  # Initialize function value (to check convergence)
  lastf <- Inf

  # For each optimization cycle...
  while (cycle <= control$ncycles) {

    # Optimize the likelihood function
    out <- this_optimizer(fun = fun, pars = pars, control = control, extra = extra)

    # Use optimized parameters as starting point for the next cycle
    pars <- out$pars
    newf <- out$fvalue

    # Compute convergence criterion
    has_conv <- abs(newf - lastf) < control$ctol

    # TODO: This tolerance is the absolute tolerance also used in downstream
    # functions in the original code.

    # Break the loop if converged
    if (has_conv) break

    # Update the latest function value
    lastf <- newf

    # Increment cycles
    cycle <- cycle + 1L

  }

  # TODO: Message saying whether convergence is good or not.

  return(out)

}
