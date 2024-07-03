## Here we test that the simplex optimization routine works.

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

  # Precise control
  control <- list(rtolx = 1e-6, rtolf = 1e-7, atolx = 1e-9)

  # Find the parameter values that maximize the function
  out <- simplex(
    biquadratic, pars = c(x = -10, y = -10),
    extra = list(a = -1, b = -1),
    control = control
  )

  # Note: with a and b being negative there is a clear maximum.

  # Check that we have found the right maximum
  expect_true(all(round(out$pars, 3L) == 1))
  expect_equal(round(out$fvalue, 6L), 2)

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
    biquadratic, pars = c(x = -10, y = -10), extra = list(a = -0.5, b = -0.5),
    control = list(delta = 0.1)
  )

  # Should have run
  expect_true(is.list(out))

})
