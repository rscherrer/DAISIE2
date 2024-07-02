## Here we test the function that says whether an object is TRUE or FALSE.

# Test the function
test_that("Is it a yes or a no?", {

  # Use cases
  expect_true(is_yes_no(TRUE))
  expect_true(is_yes_no(FALSE))

  # Should not work on anything else
  expect_false(is_yes_no(NA))
  expect_false(is_yes_no(NaN))
  expect_false(is_yes_no("Hey!"))
  expect_false(is_yes_no(3.333))

  # Not even for a vector
  expect_false(is_yes_no(c(TRUE, TRUE)))

})
