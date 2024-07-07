## Here test the likelihood function.

# Create some data
data <- list(

  blossom = make_clade(is_present = FALSE, tcol = -10),
  bubbles = make_clade(is_present = TRUE, tcol = -5, branching_times = c(-2, -1)),
  buttercup = make_clade(is_present = FALSE, tmax = -10, tmin = -5)

)

# Parameter values
pars <- list(lambda_c = 0.18, mu = 0.02, gamma = 0.02, lambda_a = 2)

# The likelihood function should work
test_that("Use case", {

  # This should work
  expect_true(is.numeric(calc_loglik(data, pars, island_age = -30, M = 100L, nmax = 10L)))

})
