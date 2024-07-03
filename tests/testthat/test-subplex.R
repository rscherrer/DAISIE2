## Here we test the subplex optimization routine.

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

# Test a regular use case
test_that("Use case", {

  out <- subplex(biquadratic, pars = c(x = -10, y = -10), extra = list(a = -1, b = -1))

  expect_true(is.list(out))
  expect_true(all(round(out$pars, 5L) == c(1, 1)))
  expect_equal(out$fvalue, 2)
  expect_equal(out$conv, 0L)

})
