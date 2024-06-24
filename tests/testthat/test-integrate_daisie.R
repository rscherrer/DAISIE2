## Here we test the integration of the partial likelihood for a given clade.

# Parameter values
pars <- list(lambda_c = 0.18, mu = 0.02, gamma = 0.02, lambda_a = 2)

# Use case with nothing on the island whatsoever
test_that("Nothing on the island", {

  # Integrate a case with no extant clade
  expect_true(is.numeric(integrate_clade(island_age = -30, pars, nmax = 10L)))

})

# Extant singleton with known colonization time
test_that("Singleton with known colonization", {

  # Integrate a case with one extant singleton
  expect_true(is.numeric(integrate_clade(island_age = -30, pars, nmax = 10L, tcol = -10)))

})

# Extant singleton with minimum colonization time
test_that("Singleton with known (minimum) colonization", {

  # Integrate a case with one extant singleton
  expect_true(is.numeric(integrate_clade(island_age = -30, pars, nmax = 10L, tmax = -20, tmin = -10)))

})

# Extant singleton with unknown (max.) colonization time
test_that("Singleton with unknown (max.) colonization", {

  # Integrate a case with one extant singleton
  expect_true(is.numeric(integrate_clade(island_age = -30, pars, nmax = 10L, tmax = -20)))

})

# Now with an established clade
test_that("Established clade", {

  # Integrate a case with an extant clade
  expect_true(is.numeric(integrate_clade(
    island_age = -30, pars, nmax = 10L, tcol = -20,
    branching_times = c(-15, -10)
  )))

})

# Now with an established clade with unknown colonization time
test_that("Established clade with unknown colonization", {

  # Integrate a case with an extant clade
  expect_true(is.numeric(integrate_clade(
    island_age = -30, pars, nmax = 10L, tmax = -20,
    branching_times = c(-15, -10)
  )))

})

# Error if the integration method does not exist
test_that("No wonky integration method", {

  # Check that
  expect_error(integrate_clade(island_age = -30, pars, nmax = 10L, method = "hello"))

})
