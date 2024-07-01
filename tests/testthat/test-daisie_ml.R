## Here we test the likelihood maximization routine.

# Create some data
data <- list(

  blossom = make_clade(is_present = FALSE, tcol = -10),
  bubbles = make_clade(is_present = TRUE, tcol = -5, branching_times = c(-2, -1)),
  buttercup = make_clade(is_present = FALSE, tmax = -10, tmin = -5)

)

# Parameter values
pars <- list(lambda_c = 0.18, mu = 0.02, gamma = 0.02, lambda_a = 2)

# Use case
test_that("Use case", {

  daisie_ml(data, pars, island_age = -30, M = 100, nmax = 10)

})
