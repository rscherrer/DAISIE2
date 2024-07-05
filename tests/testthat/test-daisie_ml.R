## Here we test the likelihood maximization routine.

# Create some data
data <- list(

  blossom = make_clade(is_present = FALSE, tcol = -10),
  bubbles = make_clade(is_present = TRUE, tcol = -5, branching_times = c(-2, -1)),
  buttercup = make_clade(is_present = FALSE, tmax = -10, tmin = -5)

)

# Parameter values
pars <- list(lambda_c = 0.18, mu = 0.02, gamma = 0.02, lambda_a = 2)

# Make sure the routine works
test_that("Maximum likelihood", {

  # Run the optimization
  out <- daisie_ml(
    data, pars, island_age = -30, M = 100L, nmax = 10L,
    control_ml = list(maxiter = 1L), method = "subplex"
  )

  # TODO: Figure where to turn island age into a positive number.

  # Check names
  expect_true(all(names(out) == c("pars", "fvalue", "conv")))

  # Should have converged
  expect_false(out$conv == 0L)

})

# Try with the simplex algorithm
test_that("With simplex", {

  # Run the optimization
  out <- daisie_ml(
    data, pars, island_age = -30, M = 100L, nmax = 10L,
    control_ml = list(maxiter = 1L), method = "simplex"
  )

  # Check names
  expect_true(all(names(out) == c("pars", "fvalue", "conv")))

  # Should have converged
  expect_false(out$conv == 0L)

})
