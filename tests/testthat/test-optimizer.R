## Here we test the optimization routine.

# Function implementing the bivariate quadratic formula
biquadratic <- function(pars, a = -1, b = -1, c = 0, d = 0, e = 0, f = 0) {

  # pars: the two variables (named x and y)
  # a...f: coefficients of the polynomial

  # Note: with negative quadratic coefficients there is a clear maximum.

  # Extract the two variables
  x <- pars["x"]
  y <- pars["y"]

  # The formula
  return(a * x^2 + b * y^2 + c * x * y + d * x + e * y + f)

}

# Example starting values
pars0 <- c(x = 10, y = 10)

# Function to re-use the tests applied to the solution
check_solution <- function(out) {

  # out: the output

  # Check that the solution is as expected
  expect_true(all(round(out$pars, 4L) == c(0, 0)))
  expect_true(round(out$fvalue, 6L) == 0)
  expect_true(out$conv == 0)

}

# Use case
test_that("Use case", {

  # Find the maximum
  out <- optimizer(biquadratic, pars0)

  # Check that the solution is as expected
  check_solution(out)

})

# TODO: Should work with unnamed parameters (that depends on what the function
# expects really).

# TODO: Also what happens when the length of parameters is not as expected
# by the function? It also sounds like a function problem.

# Error when...
test_that("Abuse cases", {

  # Function is not a function
  expect_error(optimizer("hello", pars0))

  # Parameters are not numeric
  expect_error(optimizer(biquadratic, pars = c(x = "hi", y = "hey")))
  expect_error(optimizer(biquadratic, pars = c(x = NA, y = 10)))

  # Control options are not in a named list
  expect_error(optimizer(biquadratic, pars0, control = "hello"))
  expect_error(optimizer(biquadratic, pars0, control = list(1000L)))

  # Extra parameters are not in a list
  expect_error(optimizer(biquadratic, pars0, extra = "hello"))

  # Verbose is not TRUE or FALSE
  expect_error(optimizer(biquadratic, pars0, verbose = "hello"))
  expect_error(optimizer(biquadratic, pars0, verbose = NA))

  # Wrong method
  expect_error(optimizer(biquadratic, pars0, method = 999))
  expect_error(optimizer(biquadratic, pars0, method = "magic"))

  # Number of cycles is not a positive integer
  expect_error(optimizer(biquadratic, pars0, control = list(ncycles = "hey")))
  expect_error(optimizer(biquadratic, pars0, control = list(ncycles = -1)))
  expect_error(optimizer(biquadratic, pars0, control = list(ncycles = 3.1459)))

  # Tolerance between cycles is not a positive number
  expect_error(optimizer(biquadratic, pars0, control = list(ctol = "hey")))
  expect_error(optimizer(biquadratic, pars0, control = list(ctol = -1)))

  # Absolute tolerance is not a positive number
  expect_error(optimizer(biquadratic, pars0, control = list(atol = "hey")))
  expect_error(optimizer(biquadratic, pars0, control = list(atol = -1)))

  # Relative tolerance is not a positive number
  expect_error(optimizer(biquadratic, pars0, control = list(rtol = "hey")))
  expect_error(optimizer(biquadratic, pars0, control = list(rtol = -1)))

  # Maximum number of iterations is not a positive integer
  expect_error(optimizer(biquadratic, pars0, control = list(maxiter = "hey")))
  expect_error(optimizer(biquadratic, pars0, control = list(maxiter = -1)))
  expect_error(optimizer(biquadratic, pars0, control = list(maxiter = 3.1459)))

  # Jitter is not a number
  expect_error(optimizer(biquadratic, pars0, control = list(jitter = NA)))
  expect_error(optimizer(biquadratic, pars0, control = list(jitter = "hey")))

  # Scaling factors are not numbers in a vector of the right size
  expect_error(optimizer(biquadratic, pars0, control = list(scale = NA)))
  expect_error(optimizer(biquadratic, pars0, control = list(scale = "hey")))
  expect_error(optimizer(biquadratic, pars0, control = list(scale = c(1, 1, 1, 1))))

})

# Use case with verbose
test_that("Use case with verbose", {

  # Find the maximum
  out <- suppressMessages(optimizer(biquadratic, pars0, verbose = TRUE))

  # Check that the solution is as expected
  check_solution(out)

})

# Suppress warnings
test_that("Suppress warnings", {

  # Optimize with warning suppressor
  out <- optimizer(biquadratic, pars0, warn = FALSE)

  # Check that the solution is as expected
  check_solution(out)

})

# Unknown options supplied should be ignored
test_that("Unknown options are ignored", {

  # Check that the right message is displayed
  suppressMessages(expect_message(
    optimizer(biquadratic, pars0, control = list(mojo = "jojo"), verbose = TRUE),
    regexp = "^Unknown option.* mojo\n$"
  ))

  # Same with a mix of known and unknown options
  suppressMessages(expect_message(
    optimizer(biquadratic, pars0, control = list(maxiter = 10L, mojo = "jojo"), verbose = TRUE),
    regexp = "^Unknown option.* mojo\n$"
  ))

  # Still make sure that the output is the right one
  out <- optimizer(biquadratic, pars0, control = list(mojo = "jojo"))

  # Check that the solution is as expected
  check_solution(out)

})

# Check that it works with re-scaled parameters
test_that("Re-scaled parameters", {

  # Try with a function to untransform the parameters
  out <- optimizer(biquadratic, pars = c(x = 1, y = 1), control = list(untrans = function(x) x * 10, rtol = 1e-7))

  # Check that the solution is as expected
  check_solution(out)

})

# Try multiple cycles of optimization
test_that("Multiple cycles", {

  # Still make sure that the output is the right one
  out <- suppressMessages(optimizer(biquadratic, pars0, control = list(ncycles = 10L), verbose = TRUE))

  # Note: verbose is turned on to increase coverage.

  # Check that the solution is as expected
  check_solution(out)

})

# Subplex with zero iterations
test_that("Zero iterations", {

  # Run the optimizer
  out <- optimizer(biquadratic, pars0, list(maxiter = 0L))

  # Make sure nothing was done (and no convergence)
  expect_true(all(out$pars == pars0))
  expect_false(out$conv == 0L)

})

# Problematic output
test_that("Problematic output encountered", {

  # Artificially produce a problem encountered during function evaluation
  expect_error(optimizer(function(x, y) NA, pars0))

})

# Use case with simplex
test_that("Use case with simplex", {

  # Find the maximum
  out <- optimizer(biquadratic, pars0, method = "simplex", control = list(ncycles = 10L))

  # Check that the solution is as expected
  check_solution(out)

})

# Abuse cases with simplex
test_that("Simplex abuse cases", {

  # Absolute tolerance is not a positive number
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(atol = "hey")))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(atol = -1)))

  # Relative tolerance is not a positive number
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(rtol = "hey")))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(rtol = -1)))

  # Function tolerance is not a positive number
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(ftol = "hey")))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(ftol = -1)))

  # Maximum number of iterations is not a positive integer
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(maxiter = "hey")))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(maxiter = -1)))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(maxiter = 3.1459)))

  # Parameter dodging is not a positive number
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(delta = "hey")))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(delta = -1)))

  # Dodge for zero-values is not a number
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(dzero = NA)))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(dzero = "hey")))

  # Geometric transformation parameter rho is not a positive number
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(rho = "hey")))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(rho = -1)))

  # Geometric transformation parameter chi is not a positive number
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(chi = "hey")))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(chi = -1)))

  # Geometric transformation parameter psi is not a positive number
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(psi = "hey")))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(psi = -1)))

  # Geometric transformation parameter sigma is not a positive number
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(sigma = "hey")))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(sigma = -1)))

  # Back-transformation is neither NULL nor a function (triggered from within simplex)
  expect_error(simplex(biquadratic, pars0, control = list(untrans = "hey")))
  expect_error(simplex(biquadratic, pars0, control = list(untrans = 1)))
  expect_error(simplex(biquadratic, pars0, control = list(untrans = NA)))

  # Transformation is neither NULL nor a function
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(trans = "hey")))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(trans = 1)))
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(trans = NA)))

  # Transformation functions
  trans <- function(x) x / 10
  untrans <- function(x) x * 10

  # Only one of those two functions is provided
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(trans = trans)))
  suppressMessages(expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(untrans = untrans), verbose = TRUE)))

  # Wrong back-tranformation
  wrong_untrans <- function(x) x * 5

  # Both are provided but they are no inverse of each other
  expect_error(optimizer(biquadratic, pars0, method = "simplex", control = list(trans = trans, untrans = wrong_untrans)))

})

# Simplex with verbose
test_that("Simplex with verbose", {

  # Find the maximum
  out <- suppressMessages(optimizer(biquadratic, pars0, method = "simplex", verbose = TRUE, control = list(ncycles = 10L)))

  # Check that the solution is as expected
  check_solution(out)

})

# Simplex with verbose and no convergence
test_that("Simplex with verbose and no convergence", {

  # Find the maximum
  out <- suppressMessages(optimizer(biquadratic, pars0, method = "simplex", verbose = TRUE))

  # Check that the algorithm has not converged
  expect_false(out$conv == 0L)

})

# Simplex with re-scaled parameters
test_that("Simplex with re-scaled parameters", {

  # Transformation functions
  trans <- function(x) x / 10
  untrans <- function(x) x * 10

  # Try with a function to untransform the parameters
  out <- optimizer(biquadratic, pars = c(x = 1, y = 1), method = "simplex", control = list(untrans = untrans, trans = trans, ncycles = 10L))

  # Check that the solution is as expected
  check_solution(out)

})

# Try a function that triggers the shrink operation in simplex
test_that("Trigger simplex shrinkage", {

  # Run an example that involves shrinking
  out <- optimizer(
    biquadratic,
    pars = c(x = -10, y = -10),
    method = "simplex",
    control = list(delta = 0.1),
    extra = list(a = -0.5, b = -0.5, c = 1, d = 1, e = 1, f = 1)
  )

  # Make sure it has run
  expect_true(is.list(out))

})

# TODO: Our simplex algorithm sucks.

# Simplex with one-parameter function
test_that("Simplex with one parameter", {

  # Function with one parameter
  FUN <- function(x) -x^2

  # Find its maximum
  out <- optimizer(FUN, pars = 0, method = "simplex")

  # Make sure we have found it
  expect_equal(out$pars, 0)
  expect_equal(out$fvalue, 0)
  expect_equal(out$conv, 0L)

})

# Simplex when infinity is reached
test_that("Simplex with infinity", {

  # Square function
  FUN <- function(x) x^2

  # Run the optimizer
  out <- optimizer(FUN, pars = 10, method = "simplex")

  # Hopefully the maximum is infinite
  expect_equal(out$fvalue, Inf)
  expect_equal(out$conv, 0L)

})

# TODO: Go over names of tests across the entire package.
