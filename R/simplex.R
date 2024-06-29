# Function to set the control options
make_control <- function(options, n = 1L) {

  # options: named list of user-specified options
  # n: number of parameters to optimize

  # Initial checks
  testit::assert(is.list(options))
  testit::assert(is_number(n, scalar = TRUE, integer = TRUE, sign = 1))

  # Default control parameters
  control <- list(

    # Numerical tolerance parameters to check convergence
    rtolx = 1e-4,
    rtolf = 1e-5,
    atolx = 1e-7,

    # Maximum number of iterations
    maxiter = 1000L * round((1.25)^n),

    # By how much to shift the initial guesses to initialize the simplex (e.g. 0.05 for 5%)
    delta = 0.05,

    # Value by which to shift in case initial guess is zero
    dzero = 0.00025,

    # Geometric transformation parameters
    rho = 1,
    chi = 2,
    psi = 0.5,
    sigma = 0.5

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
check_control <- function(control) {

  # control: list containing the control parameters

  # Check that it is a list
  testit::assert(is.list(control))

  # Names of the control parameters
  control_names <- c(
    "rtolx", "rtolf", "atolx", "maxiter", "delta", "dzero",
    "rho", "chi", "psi", "sigma"
  )

  # Check that they are all here
  testit::assert(all(names(control) %in% control_names))
  testit::assert(all(control_names %in% names(control)))

  # With the control parameters...
  with(control, {

    # Perform checks
    testit::assert(is_number(rtolx, sign = 1))
    testit::assert(is_number(rtolf, sign = 1))
    testit::assert(is_number(atolx, sign = 1))
    testit::assert(is_number(maxiter, integer = TRUE, sign = 1))
    testit::assert(is_number(delta, sign = 1))
    testit::assert(is_number(dzero))
    testit::assert(all(is_number(c(rho, chi, psi, sigma), scalar = FALSE, sign = 1)))

  })
}

# Wrapper around a function to avoid argument collision between environments
call_fun <- function(fun, x, extra = list()) {

  # fun: the function to call
  # x: the first argument
  # extra: list with all the other arguments

  # Call the function
  do.call(fun, c(list(x), extra))

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
get_updated_vertex <- function(x, xbar, fun, fmin, fstm, fmax, rho, chi, psi, sigma, extra = list()) {

  # x: coordinates of the worst vertex
  # xbar: centroid of the edge to flip along
  # fun: the function to measure
  # fmin: minimum (i.e. best) score out of all vertices
  # fstm: second-to-maximum score out of all vertices
  # fmax: maximum (i.e. worst) score out of all vertices
  # rho, chi, psi, sigma: update parameters
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
  fxr <- -call_fun(fun, xr, extra)

  # TODO: Can we avoid passing the function so many times?

  # TODO: Does that mean we can make it work with lists of parameters?

  # Set reflection as our new prospect
  out <- new_vertex(xr, fxr, how = "reflect")

  # If the reflection is a new minimum...
  if (fxr < fmin) {

    # Expand it to see if we can minimize further
    xe <- (1 + rho * chi) * xbar - rho * chi * x
    fxe <- -call_fun(fun, xe, extra)

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
    fxco <- -call_fun(fun, xco, extra)

    # Keep the contraction if it does better
    if (fxco <= fxr) return(new_vertex(xco, fxco, how = "contract_outside"))

    # Otherwise shrink the simplex
    return(new_vertex(NULL, NULL, how = "shrink"))

  }

  # If the reflection was the worst, then compute the inside contraction
  xci <- (1 - psi) * xbar + psi * x
  fxci <- -call_fun(fun, xci, extra)

  # Keep this contraction if it does better
  if (fxci < fmax) return(new_vertex(xci, fxci, how = "contract_inside"))

  # Otherwise shrink the simplex
  return(new_vertex(NULL, NULL, "shrink"))

}

# Function to compute (Nelder-Mead) simplex optimization
simplex <- function(

  fun, pars, control = list(), untrans = NULL, trans = NULL, extra = list()

) {

  # fun: the function (must take pars as first argument, not nec. with that name)
  # pars: vector of parameters to optimize (initial guesses, can be named)
  # control: a list of named options for the algorithm, as defined in control()
  # untrans: function to retrieve the parameters on their natural scale (the inverse of trans)
  # trans: function to return the parameters back to their transformed scale
  # extra: named list of extra arguments of the function to optimize

  # Checks
  testit::assert(is.function(fun))
  testit::assert(is_number(pars, scalar = FALSE))

  # If one transformation is provided then both must be
  testit::assert(is.function(trans) == is.function(untrans))

  # Number of parameters (i.e. dimensions)
  npars <- length(pars)

  # Update default options with user choices
  control <- make_control(control, npars)

  # Check the control parameters
  check_control(control)

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

  # Untransform if needed
  if (!is.null(untrans)) dpars <- untrans(pars)

  # Dodge
  dpars <- dpars * r0

  # Transform them back if needed
  if (!is.null(trans)) dpars <- trans(dpars)

  # Measure the transformation factors (ratio) on the transformed scale
  r <- dpars / pars

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

  # Small deviation for zero-values coordinates
  V[ii][V[ii] == 0] <- control$dzero

  # Add parameter names to the simplex if needed
  if (!is.null(names(pars))) rownames(V) <- names(pars)

  # Compute the function for each vertex
  fvalues <- apply(V, 2L, \(x) { -call_fun(fun, x, extra) })

  # TODO: Make this work with named vectors of parameters

  # Note: this algorithm finds maximizes a function by finding the minimum
  # of its negative.

  # Reorder vertices from minimum to maximum function value
  jj <- order(fvalues)
  V <- V[, jj]
  if (npars == 1L) V <- matrix(V, nrow = npars)
  fvalues <- fvalues[jj]

  # Note: this way the last vertex is the worst one.

  # For each iteration...
  for (iter in 1:control$maxiter) {

    # Compute useful metrics for convergence checking
    maxdf <- max(abs(fvalues - fvalues[1]))
    dV <- abs(V - V[, 1])

    # Note: this latter line compute the distance between each vertex (column)
    # and the first one.

    # Measure convergence in both function evaluation and simplex coordinates
    is_conf <- !is.nan(maxdf) & maxdf <= control$rtolf * abs(fvalues[1])
    is_conv <- with(control, max(dV - rtolx * abs(V[, 1])) <= 0 & max(dV) <= atolx)

    # Break the loop if convergence is reached
    if (is_conf | is_conv) break

    # Centroid of the n (not n + 1)-dimensional edge along which to flip the simplex
    xbar <- rowSums(V[, 1:npars]) / npars

    # Note: this would be a line in a 2-dimensional parameter space here
    # the simplex is a triangle, for example, but it would be a triangle
    # in a 3-dimensional parameter space where the simplex is a tetrahedron.

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

      # Extra arguments for the function to optimize
      extra

    ))

    # If the simplex must be shrunk...
    if (new_vertex$how == "shrink") {

      # For each vertex except the best one
      kk <- 2:nvertices

      # Shrink the simplex by shifting that vertex towards the best one
      V[, kk] <- V[, 1] + control$sigma * (V[, kk] - V[, 1])
      fvalues[kk] = -call_fun(fun, V[, kk], extra)

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

    # Increment iteration
    iter <- iter + 1L

  }

  # Return...
  return(list(

    # The best parameters found and their associated function value
    pars = V[, 1],
    fvalue = -fvalues[1],

    # Whether convergence has been reached
    conv = iter < control$maxiter

  ))

  # TODO: Figure why many outputs are wrapped in invisible().

  # TODO: Make simplex() its own package? Or have a precise help file and export it.

}
