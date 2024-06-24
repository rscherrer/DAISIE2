## Here test the likelihood function.

# Function to create a clade object
make_clade <- function(

  is_present = FALSE, tcol = NULL, branching_times = c(),
  tmax = NULL, tmin = NULL

) {

  # is_present: whether or not the mainland relative is on the island
  # tcol: known colonization time (could be NULL)
  # branching_times: vector of branching times (could be empty)
  # tmax: maximum colonization time (could be NULL)
  # tmin: minimum colonization time (could be NULL)

  # TODO: Might bring the checks to this function.

  # If no colonization is provided

  # Assemble into a clade object
  return(list(
    is_present = is_present,
    tcol = tcol,
    branching_times = branching_times,
    tmax = tmax,
    tmin = tmin
  ))

}

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
