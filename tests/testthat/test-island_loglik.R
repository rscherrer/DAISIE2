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
  expect_true(is.numeric(island_loglik(data, pars, island_age = -30, M = 100L, nmax = 10L)))

})

# Conditioning
test_that("Conditioning", {

  # This should work
  expect_true(is.numeric(island_loglik(data, pars, island_age = -30, M = 100L, nmax = 10L, condition = 1L)))
  expect_true(is.numeric(island_loglik(data, pars, island_age = -30, M = 100L, nmax = 10L, condition = 3L)))

})

# Abuse cases
test_that("Abuse", {

  # Should error if parameters are not positive numbers
  expect_error(island_loglik(data, pars = list(lambda_c = NA, mu = 0.02, gamma = 0.02, lambda_a = 2), island_age = -30, M = 100L, nmax = 10L))
  expect_error(island_loglik(data, pars = list(lambda_c = 0.18, mu = -0.5, gamma = 0.02, lambda_a = 2), island_age = -30, M = 100L, nmax = 10L))
  expect_error(island_loglik(data, pars = list(lambda_c = 0.18, mu = -0.5, gamma = 0.02, lambda_a = 2), island_age = -30, M = 100L, nmax = 10L, optimized = TRUE))
  expect_error(island_loglik(data, pars, island_age = -30, M = 100L, nmax = 10L, condition = 10L))

  # TODO: This last one should be superseded by a test coming directly from the
  # maximization procedure.

})

