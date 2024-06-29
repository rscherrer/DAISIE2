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
  sol <- simplex(
    biquadratic, pars = c(x = -10, y = -10),
    extra = list(a = -1, b = -1),
    control = control
  )

  # Note: with a and b being negative there is a clear maximum.

  # Check that we have found the right maximum
  expect_true(all(round(sol$pars, 3L) == 1))
  expect_equal(round(sol$fvalue, 6L), 2)

})

# TODO: Check that it should error if the names of the parameters are wrong?
