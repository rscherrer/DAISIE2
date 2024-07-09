## Here we test the likelihood maximization routine.

# Fake clades
blossom <- make_clade(is_present = FALSE, tcol = -10)
bubbles <- make_clade(is_present = TRUE, tcol = -5, branching_times = c(-2, -1))
buttercup <- make_clade(is_present = FALSE, tmax = -10, tmin = -5)

# Create some data
data <- list(blossom = blossom, bubbles = bubbles, buttercup = buttercup)

# Parameter values
pars <- list(lambda_c = 0.18, mu = 0.02, gamma = 0.02, lambda_a = 2)

# Make sure the routine works
test_that("Maximum likelihood", {

  # Run the optimization
  out <- daisie_ml(
    data, pars, island_age = -30, M = 100L, nmax = 10L,
    control_ml = list(maxiter = 0L)
  )

  # Check names
  expect_true(all(names(out) == c("pars", "fvalue", "conv")))

  # Should have converged
  expect_false(out$conv == 0L)

})

# Abuse cases
test_that("Abuse", {

  # Parameters are not a list
  expect_error(daisie_ml(data, pars = "hello", island_age = -30, M = 100L, nmax = 10L))

  # Parameters are not named
  expect_error(daisie_ml(data, pars = list(1), island_age = -30, M = 100L, nmax = 10L))

  # Unknown parameters
  expect_error(daisie_ml(data, pars = list(foo = 1, bar = 1), island_age = -30, M = 100L, nmax = 10L))

  # Missing parameters
  expect_error(daisie_ml(data, pars = list(lambda_c = 0.1), island_age = -30, M = 100L, nmax = 10L))

  # Parameters are not positive numbers
  expect_error(daisie_ml(data, list(lambda_c = -0.18, mu = 0.02, gamma = 0.02, lambda_a = 2), island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(data, list(lambda_c = "hello", mu = 0.02, gamma = 0.02, lambda_a = 2), island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(data, list(lambda_c = NA, mu = 0.02, gamma = 0.02, lambda_a = 2), island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(data, list(lambda_c = 0.18, mu = -0.02, gamma = 0.02, lambda_a = 2), island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(data, list(lambda_c = 0.18, mu = 0.02, gamma = -0.02, lambda_a = 2), island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(data, list(lambda_c = 0.18, mu = 0.02, gamma = 0.02, lambda_a = -2), island_age = -30, M = 100L, nmax = 10L))

  # Island age is not a negative number
  expect_error(daisie_ml(data, pars, island_age = 2, M = 100L, nmax = 10L))
  expect_error(daisie_ml(data, pars, island_age = "hi", M = 100L, nmax = 10L))

  # Mainland pool is not a positive integer
  expect_error(daisie_ml(data, pars, island_age = -30, M = -3L, nmax = 10L))
  expect_error(daisie_ml(data, pars, island_age = -30, M = 3.1459, nmax = 10L))

  # Possible number of unobserved species is not a positive integer
  expect_error(daisie_ml(data, pars, island_age = -30, M = 100L, nmax = -3L))
  expect_error(daisie_ml(data, pars, island_age = -30, M = 100L, nmax = 3.1459))

  # Method is not a character string
  expect_error(daisie_ml(data, pars, island_age = -30, M = 100L, nmax = 10L, method = 1))
  expect_error(daisie_ml(data, pars, island_age = -30, M = 100L, nmax = 10L, method = NA))

  # Verbose is not TRUE or FALSE
  expect_error(daisie_ml(data, pars, island_age = -30, M = 100L, nmax = 10L, verbose = NA))
  expect_error(daisie_ml(data, pars, island_age = -30, M = 100L, nmax = 10L, verbose = "hi"))

  # Control options are not in a list
  expect_error(daisie_ml(data, pars, island_age = -30, M = 100L, nmax = 10L, control_ml = "hi"))
  expect_error(daisie_ml(data, pars, island_age = -30, M = 100L, nmax = 10L, control_ml = 1))

  # Data is not a list
  expect_error(daisie_ml(data = "hi", pars, island_age = -30, M = 100L, nmax = 10L))

  # Data is empty
  expect_error(daisie_ml(data = list(), pars, island_age = -30, M = 100L, nmax = 10L))

  # Clades have no names
  expect_error(daisie_ml(data = unname(data), pars, island_age = -30, M = 100L, nmax = 10L))

  # Clade is not a list
  expect_error(daisie_ml(c(data, wrong = "hi"), pars, island_age = -30, M = 100L, nmax = 10L))

  # Clade does not have the right number of elements
  expect_error(daisie_ml(c(data, wrong = list(list())), pars, island_age = -30, M = 100L, nmax = 10L))

  # Clade does does not have named elements
  expect_error(daisie_ml(c(data, wrong = list(unname(buttercup))), pars, island_age = -30, M = 100L, nmax = 10L))

  # Clade with unknown elements
  expect_error(daisie_ml(c(data, wrong = list(list(mojo = "jojo"))), pars, island_age = -30, M = 100L, nmax = 10L))

  # Clade with missing elements
  expect_error(daisie_ml(c(data, wrong = list(list(tcol = -3))), pars, island_age = -30, M = 100L, nmax = 10L))

  # Presence is not TRUE or FALSE
  expect_error(daisie_ml(c(data, wrong = list(make_clade(is_present = "hey"))), pars, island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(c(data, wrong = list(make_clade(is_present = NA))), pars, island_age = -30, M = 100L, nmax = 10L))

  # TODO: Do we want to check for specific error messages?

  # Colonization time and upper bound cannot both be missing
  expect_error(daisie_ml(c(data, wrong = list(make_clade())), pars, island_age = -30, M = 100L, nmax = 10L))

  # Colonization time is either known or unknown
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tcol = -1, tmax = -2))), pars, island_age = -30, M = 100L, nmax = 10L))

  # An upper bound is required if there is a lower bound
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tcol = -1, tmin = -0.5))), pars, island_age = -30, M = 100L, nmax = 10L))

  # TODO: Make sure we are consistent with our use of "upper" and "lower" bound.

  # Timings are not negative numbers
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tcol = 1))), pars, island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tmax = 1))), pars, island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tmax = -1, tmin = 1))), pars, island_age = -30, M = 100L, nmax = 10L))

  # Non-negative branching times
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tcol = -1, branching_times = c(-1, 1)))), pars, island_age = -30, M = 100L, nmax = 10L))

  # Anachronisms
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tmax = -40))), pars, island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tcol = -40))), pars, island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tmax = -20, tmin = -25))), pars, island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tcol = -20, branching_times = -40))), pars, island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tcol = -20, branching_times = -25))), pars, island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tmax = -20, branching_times = -25))), pars, island_age = -30, M = 100L, nmax = 10L))
  expect_error(daisie_ml(c(data, wrong = list(make_clade(tmax = -20, tmin = -10, branching_times = -15))), pars, island_age = -30, M = 100L, nmax = 10L))

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

# TODO: Think of a case where problematic parameters are found during
# optimization (for now we very artificially produce this use case in tests).
