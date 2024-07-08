## Here we test the function we use to re-scale model parameters.

# Values to re-scale
x <- c(Inf, 0, -Inf, 0.5)

# Test the range-transformation function
test_that("Range transformation works", {

  # Constrict the range
  y1 <- range_transform(x, xmax = 1)
  y2 <- range_transform(x, xmax = 2)

  # Check values
  expect_true(all(y1 == c(1, 0, -1, 1/3)))
  expect_true(all(y2 == c(2, 0, -2, 2/3)))

  # Back-transform
  y3 <- range_transform(y1, xmax = 1, inverse = TRUE)
  y4 <- range_transform(y2, xmax = 2, inverse = TRUE)

  # Check that the back-transformation worked
  expect_true(all(round(y3, 6L) == x))
  expect_true(all(round(y4, 6L) == x))

})

# Abuse cases
test_that("Abuse", {

  # Input is not a numeric vector
  expect_error(range_transform("hey"))

  # Upper bound not a finite positive number
  expect_error(range_transform(x, xmax = NA))
  expect_error(range_transform(x, xmax = -3))
  expect_error(range_transform(x, xmax = Inf))

  # Inverse is not TRUE or FALSE
  expect_error(range_transform(x, inverse = NA))

  # Values out of range
  expect_error(range_untransform(range_transform(x), xmax = 0.5))

  # TODO: Replace the inverse of range_transform() with range_untransform()
  # wherever needed.

})

# Edge cases
test_that("Shrinking to zero", {

  # Shrink everything to zero
  expect_true(all(range_transform(x, xmax = 0) == 0))

})
