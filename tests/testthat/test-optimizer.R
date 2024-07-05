## Here we test the optimization routine.

# Function implementing the bivariate quadratic formula
biquadratic <- function(pars, a = 1, b = 1, c = 1, d = 1, e = 1, f = 1) {

  # pars: the two variables (named x and y)
  # a...f: coefficients of the polynomial

  # Extract the two variables
  x <- pars["x"]
  y <- pars["y"]

  # The formula
  return(a * x^2 + b * y^2 + c * x * y + d * x + e * y + f)

}

# Check that we can find the maximum of a well-known function
test_that("Simplex on bivariate quadratic formula", {

  # Find the parameter values that maximize the function
  out <- simplex(
    biquadratic,
    pars = c(x = -10, y = -10),
    extra = list(a = -1, b = -1),
    control = list(rtol = 1e-6, ftol = 1e-7, atol = 1e-9)
  )

  # Note: with a and b being negative there is a clear maximum.

  # Check that we have found the right maximum
  expect_true(all(round(out$pars, 3L) == 1))
  expect_equal(round(out$fvalue, 6L), 2)

})

# Check that parameter re-scaling works
test_that("Parameter re-scaling", {

  # Find the parameter values that maximize the function
  out <- simplex(
    biquadratic,
    pars = c(x = -10, y = -10) / 10,
    extra = list(a = -1, b = -1),
    control = list(trans = \(x) x / 10, untrans = \(x) x * 10, rtol = 1e-6, ftol = 1e-7, atol = 1e-9)
  )

  # Check that we have found the right maximum
  expect_true(all(round(10 * out$pars, 3L) == 1))

})

# Check that parameter re-scaling works
test_that("Subplex with parameter re-scaling", {

  # Find the parameter values that maximize the function
  out <- subplex(
    biquadratic,
    pars = c(x = -10, y = -10) / 10,
    extra = list(a = -1, b = -1),
    control = list(untrans = \(x) x * 10)
  )

  # Check that we have found the right maximum
  expect_true(all(round(10 * out$pars, 3L) == 1))

})

# Error if wrong parameters
test_that("Wrong parameters", {

  # Should error
  expect_error(simplex(biquadratic, pars = c(hey = 2, hello = 2)))

})

# Can handle infinite optima
test_that("Infinite optimum", {

  # Find the solution of a concave polynomial
  out <- simplex(biquadratic, pars = c(x = -10, y = -10))

  # Check that the maximum is infinite
  expect_equal(out$fvalue, Inf)

})

# Function with only one parameter
test_that("Function with one parameter", {

  # Find the solution
  out <- simplex(fun = \(x) -x^2, pars = c(x = 0))

  # Make sure the maximum was found
  expect_equal(out$pars, 0)
  expect_equal(out$fvalue, 0)

})

# Test an example where the simplex shrinks
test_that("Shrinking simplex", {

  # Run an example where the simplex has to shrink at some point
  out <- simplex(
    biquadratic,
    pars = c(x = -10, y = -10),
    extra = list(a = -0.5, b = -0.5),
    control = list(delta = 0.1)
  )

  # Should have run
  expect_true(is.list(out))

})

# Test a regular use case
test_that("Use case", {

  # Optimize
  out <- subplex(
    biquadratic,
    pars = c(x = -10, y = -10),
    extra = list(a = -1, b = -1)
  )

  # Check
  expect_true(is.list(out))
  expect_true(all(round(out$pars, 5L) == c(1, 1)))
  expect_equal(out$fvalue, 2)
  expect_equal(out$conv, 0L)

})

# Check that we can find the maximum of a well-known function
test_that("Optimizer on bivariate quadratic formula", {

  # Find the parameter values that maximize the function
  out <- optimizer(
    biquadratic,
    pars = c(x = -10, y = -10),
    extra = list(a = -1, b = -1),
    method = "simplex",
    control = list(rtol = 1e-6, ftol = 1e-7, atol = 1e-9, ncycles = 100L)
  )

  # Check that we have found the right maximum
  expect_true(all(round(out$pars, 3L) == 1))
  expect_equal(round(out$fvalue, 6L), 2)

})

# Check control option setters
test_that("Control options", {

  # Example option
  options <- list(maxiter = 666L)

  # For simplex or subplex
  control1 <- make_control(options, method = "simplex")
  control2 <- make_control(options, method = "subplex")

  # Check that the right options are generated
  expect_true(all(names(control1) == names(get_default_options_simplex())))
  expect_true(all(names(control2) == names(get_default_options_subplex())))

  # For the optimizer as well
  control3 <- make_control(options, method = "simplex", meta = TRUE)

  # Check that the options contain the optimizer-specific ones
  expect_true(all(names(get_default_options_optimizer()) %in% names(control3)))

  # Check values have been correctly updated
  expect_true(all(c(control1$maxiter, control2$maxiter, control3$maxiter) == 666L))

  # Try with unknown option
  control4 <- make_control(list(mojo = "jojo"), method = "simplex")

  # Make sure it was not added
  expect_false("mojo" %in% names(control4))

})
