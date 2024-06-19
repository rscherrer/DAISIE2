## Here we check the data preparation routine.

# Example data
data <- data.frame(

  Clade_name = c("foo", "bar"),
  Status = c("Endemic", "NonEndemic"),
  Missing_species = c(3, 0),
  Branching_times = c("3, 2, 1", NA)

)

# Check that the main preparation function works
test_that("Use case", {

  # Create a DAISIE list
  datalist <- DAISIE_dataprep(data, island_age = 4, M = 100)

  # Check the first element of the list
  expect_length(datalist, 3L)
  expect_equal(datalist[[1]]$island_age, 4)
  expect_equal(datalist[[1]]$not_present, 98L)

  # Check the data for the first colonist clade
  expect_equal(datalist[[2]]$colonist_name, "foo")
  expect_true(all(datalist[[2]]$branching_times == c(4, 3, 2, 1)))
  expect_equal(datalist[[2]]$stac, 2L)
  expect_equal(datalist[[2]]$missing_species, 3L)
  expect_equal(datalist[[2]]$type1or2, 1L)

  # Check the data for the second colonist clade
  expect_equal(datalist[[3]]$colonist_name, "bar")
  expect_true(all(datalist[[3]]$branching_times == c(4, 3.99999)))
  expect_equal(datalist[[3]]$stac, 1L)
  expect_equal(datalist[[3]]$missing_species, 0L)
  expect_equal(datalist[[3]]$type1or2, 1L)

})

# Check with an input that is not a data frame
test_that("Input is not a data frame", {

  # Should error
  expect_error(DAISIE_dataprep("hello", island_age = 4, M = 100))

})

# Should error if input data has no rows
test_that("Input has no rows", {

  # Should error
  expect_error(DAISIE_dataprep(data.frame(), island_age = 4, M = 100))

})

# Should error if missing columns
test_that("Missing columns", {

  # Example data
  data <- data.frame(Clade_name = c("foo", "bar"))

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Should error if clade names are not character
test_that("Clade name not character", {

  # Example data
  data$Clade_name <- c(1, 2)

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Should error if duplicate clade names
test_that("Duplicate clade names", {

  # Example data
  data$Clade_name <- c("foo", "foo")

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Should error if unknown clade statuses
test_that("Unknown statuses", {

  # Example data
  data$Status <- c("Endemic", "Bob")

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Should error if non-numeric missing species
test_that("Non-numeric missing species", {

  # Example data
  data$Missing_species <- c("hello", 0)

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Should error if non-integral missing species
test_that("Non-integral missing species", {

  # Example data
  data$Missing_species <- c(1.1, 0)

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Should error if non-natural missing species
test_that("Non-natural missing species", {

  # Example data
  data$Missing_species <- c(-1L, 0L)

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Should work if only NAs are provided as branching times
test_that("NA-only branching times", {

  # Example data
  data$Branching_times <- c(NA, NA)

  # Should work
  datalist <- DAISIE_dataprep(data, island_age = 4, M = 100)

  # Check that colonisation times have been inferred and capped
  expect_equal(datalist[[2]]$branching_times[2], 3.99999)
  expect_equal(datalist[[3]]$branching_times[2], 3.99999)

})

# Should error if branching times are not characters or NA
test_that("Non-character branching times", {

  # Example data
  data$Branching_times <- c(42, NA)

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Should error if some branching times are unreadable
test_that("Unreadable branching times", {

  # Introduce some unexpected branching time
  data$Branching_times[1] <- "3, 2, NA, 1"

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Should error if some branching times are negative
test_that("Negative branching times", {

  # Introduce some unexpected branching time
  data$Branching_times[1] <- "3, 2, -1, 1"

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Use case with colonisation time older than the island, capped to island age
test_that("Colonisation older than island age", {

  # Introduce colonisation time older than the island
  data$Branching_times[1] <- "5, 2, 1"

  # Should work
  datalist <- DAISIE_dataprep(data, island_age = 4, M = 100)

  # Colonisation time should have been capped
  expect_equal(datalist[[2]]$branching_times[2], 3.99999)

})

# Should error if branching time(s) older than island age
test_that("Cladogenesis events older than colonisation or island age", {

  # Introduce branching time older than the island
  data$Branching_times[1] <- "6, 5, 2, 1"

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Should error if island age is not of length one
test_that("Non-scalar island age", {

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = c(1, 1), M = 100))

})

# Should error if island age is not a number
test_that("Non-numeric island age", {

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = "hello", M = 100))

})

# Should error if island age is negative
test_that("Negative island age", {

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = -4, M = 100))

})

# Should error if mainland pool is not of length one
test_that("Non-scalar mainland pool", {

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = c(1, 1)))

})

# Should error if mainland pool is not numeric
test_that("Non-numeric mainland pool", {

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = "hello"))

})

# Should error if mainland pool is not integral
test_that("Non-integral mainland pool", {

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 1.3))

})

# Should error if mainland pool is not natural
test_that("Non-natural mainland pool", {

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = -10L))

})

# Test use case with type 2 species
test_that("With type 2 species", {

  # Prepare data
  datalist <- DAISIE_dataprep(data, island_age = 4, M = 100, list_type2_clades = "foo")

  # Check attributes relevant to type 2 species
  expect_equal(datalist[[1]]$not_present_type1, 49L)
  expect_equal(datalist[[1]]$not_present_type2, 49L)
  expect_equal(datalist[[2]]$type1or2, 2L)
  expect_equal(datalist[[3]]$type1or2, 1L)

})

# Should error if unknown type 2 clades
test_that("Unknown type 2 clades", {

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100, list_type2_clades = "hello"))

})

# Test use case with custom type 2 proportion on the mainland
test_that("Custom type 2 proportion", {

  # Prepare data
  datalist <- DAISIE_dataprep(
    data, island_age = 4, M = 100,
    list_type2_clades = "foo",
    prop_type2_pool = 0.2
  )

  # Check that the custom proportion is carried over into the output
  expect_equal(datalist[[1]]$not_present_type1, 79L)
  expect_equal(datalist[[1]]$not_present_type2, 19L)

})

# Should error if type 2 proportion is not a scalar
test_that("Non-scalar type 2 proportion", {

  # Should error
  expect_error(
    DAISIE_dataprep(
      data, island_age = 4, M = 100,
      list_type2_clades = "foo",
      prop_type2_pool = c(1, 1)
    )
  )

})

# Should error if type 2 proportion is not numeric
test_that("Non-numeric type 2 proportion", {

  # Should error
  expect_error(
    DAISIE_dataprep(
      data, island_age = 4, M = 100,
      list_type2_clades = "foo",
      prop_type2_pool = "hello"
    )
  )

})

# Should error if type 2 proportion is smaller than zero
test_that("Non-numeric type 2 proportion", {

  # Should error
  expect_error(
    DAISIE_dataprep(
      data, island_age = 4, M = 100,
      list_type2_clades = "foo",
      prop_type2_pool = -1
    )
  )

})

# Should error if type 2 proportion is larger than one
test_that("Non-numeric type 2 proportion", {

  # Should error
  expect_error(
    DAISIE_dataprep(
      data, island_age = 4, M = 100,
      list_type2_clades = "foo",
      prop_type2_pool = 2
    )
  )

})

# Should error if too many branching times for non-endemic clade
test_that("Too many branching times for non-endemic", {

  # Add many branching times to the non-endemic clade
  data$Branching_times[2] <- "3, 2, 1"

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})

# Should error if non-endemic clade has missing species
test_that("Non-endemic clade with missing species", {

  # Add missing species to the non-endemic clade
  data$Missing_species[2] <- 2

  # Should error
  expect_error(DAISIE_dataprep(data, island_age = 4, M = 100))

})
