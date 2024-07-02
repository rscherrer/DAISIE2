## Tests of the functionalities of a utility function

# Test the function
test_that("Is it a number and what kind?", {

  # Check that zero is a number
  expect_true(is_number(0))

  # Check that a word is not a number
  expect_false(is_number("hello"))

  # Check that a vector of zeros is a vector of numbers
  expect_true(all(is_number(c(0, 0), scalar = FALSE)))

  # Check that a vector of numbers is not a scalar
  expect_false(is_number(c(0, 0), scalar = TRUE))

  # Check that one is an integer
  expect_true(is_number(1, integer = TRUE))

  # Check that one half is not
  expect_false(is_number(0.5, integer = TRUE))

  # Check that one is positive and not negative
  expect_true(is_number(1, sign = 1))
  expect_false(is_number(1, sign = -1))

  # Check that minus one is negative and not positive
  expect_true(is_number(-1, sign = -1))
  expect_false(is_number(-1, sign = 1))

  # Check that a nonzero number is indeed nonzero
  expect_false(is_number(1, sign = 0))
  expect_false(is_number(-1, sign = 0))

  # Only zero is not nonzero
  expect_true(is_number(0, sign = 0))

  # Check that zero is both positive and negative
  expect_true(is_number(0, sign = -1))
  expect_true(is_number(0, sign = 1))

  # But not strictly so
  expect_false(is_number(0, sign = -1, strict = TRUE))
  expect_false(is_number(0, sign = 1, strict = TRUE))

  # Check with a custom vector of numbers
  expect_true(all(is_number(c(0, -3, 2, -4, 5), sign = c(1, -1, 1, -1, 1), scalar = FALSE)))
  expect_false(all(is_number(c(0, -3, 2, -4, -5), sign = c(1, -1, 1, -1, 1), scalar = FALSE)))

  # Check variations of that function
  expect_true(is_positive(3))
  expect_false(is_positive(-3))
  expect_true(is_negative(-3))
  expect_false(is_negative(3))

  # Same for integers
  expect_true(is_positive_integer(3))
  expect_false(is_positive_integer(-3))
  expect_false(is_positive_integer(3.3))

  # And vectors
  expect_true(all(is_positive_vector(c(3, 3))))
  expect_false(all(is_positive_vector(c(3, -3))))
  expect_false(all(is_negative_vector(c(3, -3))))
  expect_true(all(is_negative_vector(-3, -3)))

})
