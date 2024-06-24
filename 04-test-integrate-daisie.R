## Here we test the main integration function.

# Parameter values
pars <- list(lambda_c = 0.18, mu = 0.02, gamma = 0.02, lambda_a = 2)

# Use case with nothing on the island whatsoever
test_that("Nothing on the island", {

  # Integrate a case with no extant clade
  out <- integrate_daisie(island_age = -30, pars, nmax = 10L)

  # Check that we end with one extant singleton
  expect_equal(out$k, 0L)

})

# Now with one extant singleton and known colonization time
test_that("Known colonization", {

  # Integrate a case with one extant singleton
  out <- integrate_daisie(island_age = -30, pars, nmax = 10L, tcol = -10)

  # Check that we end with one extant singleton
  expect_equal(out$k, 1L)

})

# Now with unknown colonization time
test_that("Unknown colonization", {

  # Integrate a case with one extant singleton
  out <- integrate_daisie(island_age = -30, pars, nmax = 10L, tmax = -20, tmin = -1)

  # Check that we end with one extant singleton
  expect_equal(out$k, 1L)

})

# Now with an established clade
test_that("Established clade", {

  # Integrate a case with an extant clade
  out <- integrate_daisie(
    island_age = -30, pars, nmax = 10L, tcol = -20,
    branching_times = c(-15, -10)
  )

  # Check that we end with one extant singleton
  expect_equal(out$k, 3L)

})

# Now an establishing clade with unknown (maximum) colonization time
test_that("Established clade with unknown colonization", {

  # Note: this makes the dynamics transitioning from k = 0 to k = 2 in one
  # go, which is a specific path in the code not covered by previous tests.

  # Integrate a case with an extant clade
  out <- integrate_daisie(
    island_age = -30, pars, nmax = 10L, tmax = -20,
    branching_times = c(-15, -10)
  )

  # Check that we end with one extant singleton
  expect_equal(out$k, 3L)

})
