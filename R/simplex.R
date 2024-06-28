# Function to update the worst vertex in the simplex
update_vertex <- function(x, xbar, fun, fmin, fstm, fmax, rho, chi, psi, sigma) {

  # x: coordinates of the worst vertex
  # xbar: centroid of the edge to flip along
  # fun: the function to measure
  # fmin: minimum (i.e. best) score out of all vertices
  # fstm: second-to-maximum score out of all vertices
  # fmax: maximum (i.e. worst) score out of all vertices
  # rho, chi, psi, sigma: update parameters

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
  fxr <- -fun(xr, ...)

  # Set reflection as our new prospect
  out <- list(x = xr, fv = fxr, how = "reflect")

  # If the reflection is a new minimum...
  if (fxr < fmin) {

    # Expand it to see if we can minimize further
    xe <- (1 + rho * chi) * xbar - rho * chi * x
    fxe <- -fun(xe, ...)

    # Return the expansion if smaller
    if (fxe < fxr) return(list(x = xe, fv = fxe, how = "expand"))

    # Otherwise stick with the reflection
    return(out)

  }

  # Else, stick with the reflection if it is better than the worst two
  if (fxr < fstm) return(out)

  # Else, if it is at least not the worst...
  if (fxr < fmax) {

    # Contract the worst vertex
    xco <- (1 + psi * rho) * xbar - psi * rho * x
    fxco <- -fun(xco, ...)

    # Keep the contraction if it does better
    if (fxco <= fxr) return(list(x = xco, fv = fxco, how = "contract_outside"))

    # Otherwise shrink the simplex
    return(list(x = NULL, fv = NULL, how = "shrink"))

  }

  # If the reflection was the worst, then compute the inside contraction
  xci <- (1 - psi) * xbar + psi * x
  fxci <- -fun(xci, ...)

  # Keep this contraction if it does better
  if (fxci < fmax) return(list(x = xci, fv = fxci, how = "contract_inside"))

  # Otherwise shrink the simplex
  return(list(x = NULL, fv = NULL, how = "shrink"))

}

# Function to compute (Nelder-Mead) simplex optimization
simplex <- function(

  fun, pars, rtolx, rtolf, atolx, maxiter, delta = 0.05, dzero = 0.00025,
  rho = 1, chi = 2, psi = 0.5, sigma = 0.5, untrans = NULL, trans = NULL, ...

) {

  # fun: the function (must take pars as first argument, not nec. with that name)
  # pars: vector of parameters to optimize
  # rtolx, rtolf, atolx: numerical tolerance parameters
  # maxiter: maximum number of iterations
  # delta: the fraction by which to shift parameters to initialize the simplex (e.g. 0.05 for 5%)
  # dzero: the value by which to dodge parameter values that are equal to zero
  # rho, chi, psi, sigma: parameters of the geometrical transformation of the simplex
  # untrans: function to retrieve the parameters on their natural scale (the inverse of trans)
  # trans: function to return the parameters back to their transformed scale
  # ...: extra arguments of the function to optimize

  # TODO: Replace sign with positive in use cases of is_number.

  # Checks
  testit::assert(is.function(fun))
  testit::assert(is_number(pars, scalar = FALSE))
  testit::assert(is_number(rtolx, sign = 1))
  testit::assert(is_number(rtolxf, sign = 1))
  testit::assert(is_number(atolx, sign = 1))
  testit::assert(is_number(maxiter, integer = TRUE, sign = 1))
  testit::assert(is_number(delta, sign = 1))
  testit::assert(is_number(dzero))
  testit::assert(all(is_number(c(rho, chi, psi, sigma), scalar = FALSE, sign = 1)))

  # If one transformation is provided both must be
  testit::assert(is.function(trans) == is.function(untrans))

  # Number of parameters (i.e. dimensions)
  n <- length(pars)

  # Number of vertices
  nvertices <- n + 1L

  # Prepare the simplex as a matrix of coordinates for each vertex
  v <- matrix(rep(pars, nvertices), nrow = n)

  # Ratio by which to dodge each vertex dimension
  r0 <- 1 + delta

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
  ii <- seq(n)
  ii <- ii * n + ii

  # Note: those indices are positions in the matrix corresponding to (i, i + 1)
  # for as many i's as there are dimensions (parameters). This way, each
  # vertex is shifted by a small amount in only one particular direction.
  # See https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method for
  # a description of the algorithm.

  # Shift each vertex by the specified ratio along only one dimension
  v[ii] <- r * v[ii]

  # Small deviation for zero-values coordinates
  v[ii][v[ii] == 0] <- dzero

  # Compute the function for each vertex
  fv <- apply(v, 2L, \(x) { -fun(x, ...) })

  # Note: this algorithm finds maximizes a function by finding the minimum
  # of its negative.

  # Reorder vertices from minimum to maximum function value
  jj <- order(fv)
  v <- v[, jj]
  if (n == 1L) v <- matrix(v, nrow = n)
  fv <- fv[jj]

  # Note: this way the last vertex is the worst one.

  # For each iteration...
  for (iter in 1:maxiter) {

    # Compute useful metrics for convergence checking
    maxdf <- max(abs(fv - fv[1]))
    dv <- abs(v - v[, 1])

    # Note: this latter line compute the distance between each vertex (column)
    # and the first one.

    # Measure convergence in both function evaluation and simplex coordinates
    is_conf <- !is.nan(maxdf) & maxdf <= rtolf * abs(fv[1])
    is_conv <- max(dv - rtolx * abs(v2)) <= 0 & max(dv) <= atolx

    # Break the loop if convergence is reached
    if (is_conf | is_conv) break

    # Centroid of the n (not n + 1)-dimensional edge along which to flip the simplex
    xbar <- rowSums(v[, 1:n]) / n

    # Note: this would be a line in a 2-dimensional parameter space here
    # the simplex is a triangle, for example, but it would be a triangle
    # in a 3-dimensional parameter space where the simplex is a tetrahedron.

    # Figure how to update the worst vertex
    newv <- update_vertex(
      v[, nvertices], xbar, fun,
      fmin = fv[1], fstm = fv[n], fmax = fv[np1],
      rho, chi, psi, sigma
    )

    # If the simplex must be shrunk...
    if (newv$how == "shrink") {

      # For each vertex except the best one
      kk <- 2:nvertices

      # Shrink the simplex by shifting that vertex towards the best one
      v[, kk] <- v[, 1] + sigma * (v[, kk] - v[, 1])
      fv[kk] = -fun(pars = v[, kk], ...)

    } else {

      # Otherwise update the worst vertex
      v[, nvertices] <- newv$x
      fv[nvertices] <- newv$fv

    }

    # Reorder vertices from minimum to maximum function value
    jj <- order(fv)
    v <- v[, jj]
    if (n == 1L) v <- matrix(v, nrow = n)
    fv <- fv[jj]

    # Increment iteration
    iter <- iter + 1L

  }

  # Return...
  return(list(

    # The best parameters found and their associated function value
    pars = v[, 1],
    fvalue = -fv[1],

    # Whether convergence has been reached
    conv = iter < maxiter

  ))

  # TODO: Figure why many outputs are wrapped in invisible().

  # TODO: Make simplex() its own package? Or have a precise help file and export it.

}
