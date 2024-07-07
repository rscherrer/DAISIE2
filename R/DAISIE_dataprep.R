## This script contains the main data preparation function.

##### Accessory functions #####

# Note: Those are used to check arguments mainly, and throw error messages
# when necessary.

# Function to tell if a numeric (vector) is integral
is_integer <- function(x, natural = FALSE) {

  # x: (hopefully) a number or a vector of numbers
  # natural: whether to only count positive integers

  # Note: it is not as simple as using is.integer() because we want e.g. 1.0
  # to be recognized as integer here.

  # Early exit if non-numeric
  if (!is.numeric(x)) return(FALSE)

  # Early exit if negative
  if (natural & any(x < 0)) return(FALSE)

  # x is an integer if its ceiling (or its floor) is itself
  all(ceiling(x) == x)

}

# Function to error with a list of missing or unknown elements
stop_wrong_content <- function(weirdos, prelude = "missing elements") {

  # weirdos: elements that we must tell the user are missing/unknown
  # prelude: beginning of the error message

  # Error with the right message
  error_msg <- paste(weirdos, collapse = ", ")
  error_msg <- paste0(prelude, ": ", error_msg)
  stop(error_msg)

}

# Function to check the column with clade names
check_clade_names <- function(clade_names) {

  # clade_names: the clade names

  # Allow character and numeric
  if (!("character" %in% class(clade_names)))
    stop("column Clade_name must be of type character")

  # Look for duplicate clade names
  counts <- table(clade_names)
  counts <- counts[counts > 1]

  # There should not be any
  if (length(counts) > 0)
    stop_wrong_content(names(counts), "duplicated clade name(s):")

}

# Function to check the column with clade statuses
check_statuses <- function(statuses) {

  # statuses: the statuses

  # Authorized status labels
  expected_statuses <- c(
    "NonEndemic", "NonEndemic_MaxAge", "Endemic", "Endemic_MaxAge",
    "Endemic_MaxAge_MinAge", "NonEndemic_MaxAge_MinAge",
    "Endemic_and_NonEndemic"
  )

  # Supplied
  statuses <- unique(statuses)

  # Any unknown status?
  wrong_statuses <- statuses[!(statuses %in% expected_statuses)]

  # If so, error
  if (length(wrong_statuses) > 0)
    stop_wrong_content(wrong_statuses, "unknown status(es) in column Status:")

}

# Function to check the column with numbers of missing species for each clade
check_missing_species <- function(missing_species) {

  # missing_species: the numbers of missing species per clade

  # Make sure that the number of missing species is a non-negative integer
  if (!is_integer(missing_species, natural = TRUE))
    stop("column Missing_species must contain only positive integers")

}

# Function to check the column with branching times
check_branching_times <- function(branching_times) {

  # branching_times: the branching times (in string form)

  # NAs allowed
  if (all(is.na(branching_times))) return()

  # Check that branching times are character strings
  if (class(branching_times) != "character")
    stop("branching times must be provided as character strings")

}

# Function to check the main input table
check_datatable <- function(data) {

  # data: the main input table

  # Check that the input is a table
  if (!("data.frame" %in% class(data)))
    stop("datatable must be of class data.frame")

  # Make sure there are clades
  if (nrow(data) == 0) stop("please provide at least one row")

  # Expected columns in the input table
  col_names <- c("Clade_name", "Status", "Missing_species", "Branching_times")

  # Are there any missing?
  missing_cols <- col_names[!(col_names %in% colnames(data))]

  # If so, error
  if (length(missing_cols) > 0)
    stop_wrong_content(missing_cols, "missing column(s) in datatable:")

  # Check each column
  check_clade_names(data[["Clade_name"]])
  check_statuses(data[["Status"]])
  check_missing_species(data[["Missing_species"]])
  check_branching_times(data[["Branching_times"]])

}

# Function to check the island age
check_island_age <- function(island_age) {

  # island_age: the island age

  # Prepare the error message
  error_msg <- "island_age must be a positive number"

  # Make sure island age is a positive number
  if (length(island_age) != 1L) stop(error_msg)
  if (!is.numeric(island_age)) stop(error_msg)
  if (island_age < 0) stop(error_msg)

}

# Function to check the size of the mainland pool
check_M <- function(M) {

  # M: the mainland pool

  # Prepare the error message
  error_msg <- "M must be a positive integer"

  # Make sure the mainland pool is a positive integer
  if (length(M) != 1L) stop(error_msg)
  if (!is_integer(M, natural = TRUE)) stop(error_msg)

}

# Check the list of type 2 clades
check_list_type2_clades <- function(list_type2_clades, clades) {

  # list_type2_clades: the list
  # clades: the known clades in the data list

  # Early exit if empty
  if (length(list_type2_clades) == 0) return()

  # Just in case there are duplicates (allowed)
  list_type2_clades <- unique(list_type2_clades)

  # Find unknown clades among type 2 clades provided
  unknown_clades <- list_type2_clades[!(list_type2_clades %in% clades)]

  # Error if any
  if (length(unknown_clades) > 0)
    stop_wrong_content(unknown_clades, "unknown type 2 clade(s):")

}

# Function to check the proportion of type 2 clades on the mainland
check_prop_type2_pool <- function(prop_type2_pool) {

  # prop_type2_pool: the proportion

  # Early exit if not provided
  if (is.null(prop_type2_pool)) return()

  # Prepare error message
  error_msg <- "prop_type2_pool must be a number between 0 and 1"

  # Check that it is a number between 0 and 1
  if (!is.numeric(prop_type2_pool)) stop(error_msg)
  if (length(prop_type2_pool) != 1L) stop(error_msg)
  if (prop_type2_pool > 1 | prop_type2_pool < 0) stop(error_msg)

}

# Function to say that unknown colonisation time will be reset to island age
message_unknown_col_time <- function(clade) {

  # clade: name of the clade

  # Display the message
  message(paste0(

    "Unknown colonisation time for ", clade,
    ". Will be set to island age."

  ))
}

# Function to say that colonisation time is older than the island
message_col_time_too_old <- function(time, clade) {

  # time: the time that is too old
  # clade: for which clade?

  # Display the message
  message(paste(

    "Colonisation time of", time, "for", clade,
    "is older than island age. Will be set to island age."

  ))

}

# Function to error because too old branching times were supplied
error_too_old_brts <- function(clade) {

  # clade: name of the clade

  # Error
  stop(paste(

    "cladogenetic event in", clade,
    "is older than the island, or of the same age as the island"

  ))

}

# Function to extract branching times for a given clade
get_branching_times <- function(

    branching_times, island_age, epss = 1e-5, verbose = FALSE, clade = ""

) {

  # branching_times: a character string encoding the branching times
  # island_age: the island age
  # epss: small offset to island age
  # verbose: whether to display messages
  # clade: clade name (for messages)

  # NA is allowed
  if (is.na(branching_times)) branching_times <- "NA"

  # If not NA, by now branching times should be in character form
  assertthat::assert_that(is.character(branching_times))

  # Separate the branching times
  branching_times <- unlist(strsplit(branching_times, split = ","))

  # Convert colonization time to NA if needed
  if (branching_times[1] == "NA") branching_times[1] <- NA

  # Count the number of NAs (maximum 1)
  n_nas <- as.integer(is.na(branching_times[1]))

  # Convert to numbers
  branching_times <- suppressWarnings(as.numeric(branching_times))

  # Error if some branching times were converted to NAs (i.e. were not numbers)
  if (sum(is.na(branching_times)) > n_nas)
    stop("some branching time(s) could not be converted to numeric")

  # Make sure all branching times are positive
  if (any(branching_times[!is.na(branching_times)] < 0))
    stop("branching times must be positive (amount of time ago)")

  # Sort the branching times in decreasing order with any NA at the front
  branching_times <- sort(branching_times, decreasing = TRUE, na.last = FALSE)

  # If the colonisation time is unknown...
  if (is.na(branching_times[1])) {

    # Set it to island age
    branching_times[1] <- island_age

    # And if needed, say it
    if (verbose) message_unknown_col_time(clade)

  }

  # If the colonisation time is older than the island...
  if (branching_times[1] > island_age & verbose) {

    # Say it
    message_col_time_too_old(branching_times[1], clade)

  }

  # Cap the colonisation time to just below island age
  branching_times[1] <- min(branching_times[1], island_age - epss)

  # By now branching times should be in strict decreasing order
  if (any(diff(branching_times) >= 0)) error_too_old_brts(clade)

  # Note: if not, some branching events are older than
  # the maximum colonisation time, and that is not allowed.

  ## TODO: can there be duplicates in branching times?

  # Append island age to the branching times
  branching_times <- c(island_age, branching_times)

  # There should by now be at least two elements (island age and col. time)
  assertthat::assert_that(length(branching_times) > 1)

  return(branching_times)

}

# Function to error because too many branchings were supplied
error_too_many_brts_for_nonendemic <- function(clade) {

  # clade: name of the clade

  # Error
  stop(paste(

    "only one branching time should be provided for", clade,
    "because it is a non-endemic species. If you mean to specifiy a",
    "minimum age as well, please use NonEndemic_MaxAge_MinAge."

  ))

}

# Function to error because too many missing species
error_too_many_missing_species_for_nonendemic <- function(clade) {

  # Error
  stop(paste(

    "Missing species for", clade, "should be 0 because it is a",
    "non-endemic species"

  ))

}

##### Main function #####

#' Prepare colonisation and branching time data to run in DAISIE
#'
#' This function produces a data object that can be run in DAISIE likelihood
#' computation/optimization functions. The function converts a user-specified
#' table to a DAISIE-compatible format. (See Galapagos_datatable.Rdata
#' for a template of an input table.)
#'
#' The output is an R list containing the data formatted to be run on other
#' DAISIE functions.
#'
#' @param datatable Data frame (table) with user-specified data. See file
#' Galapagos_datatable.Rdata for a template of an input table. Each row on the
#' table represents and independent colonisation event. The table has the following
#' four columns:
#'
#' \describe{
#'
#'  \item{\code{Clade_name}}{Name of independent colonisation event (as character)}
#'
#'  \item{\code{Status}}{
#'
#'    One of the following categories:
#'
#'    \describe{
#'
#'      \item{\code{"NonEndemic"}}{applies to non-endemic species when an approximate colonisation time is known}
#'
#'      \item{\code{"NonEndemic_MaxAge"}}{applies to non-endemic species for cases where
#'      colonisation time is unknown}
#'
#'      \item{\code{"Endemic"}}{applies to endemic species or endemic clades when an
#'      approximate colonisation time is known}
#'
#'      \item{\code{"Endemic_MaxAge"}}{applies to endemic
#'      species or endemic clades for cases where the colonisation time is unknown, or when
#'      the user wants to specify an upper bound for colonisation.
#'      This could for example apply to endemic species that have recently gone extinct because
#'      of anthropogenic causes, and which are not included
#'      in the phylogeny (\code{NA} should be given in the branching times column). It
#'      could also apply to insular radiations with long stem branches, for which the
#'      time of the first cladogenetic event is known, but the precise time of colonisation
#'      is not.}
#'
#'      \item{\code{"Endemic_MaxAge_MinAge"}}{same as \code{"Endemic_MaxAge"} but also includes a minimum
#'      age for colonisation.}
#'
#'      \item{\code{"NonEndemic_MaxAge_MinAge"}}{same as \code{"NonEndemic_MaxAge"} but
#'      also includes a minimum age for colonisation}
#'
#'      \item{\code{"Endemic_and_NonEndemic"}}{when endemic
#'      clade is present and its mainland ancestor has re-colonized}
#'
#'    }
#'  }
#'
#'  \item{\code{Missing_species}}{Number of island species that were not sampled for that
#'  particular clade (only applicable for \code{"Endemic"} clades). If \code{NA} is given in the
#'  branching times column, this should be equal to the number of species in the
#'  clade minus 1.}
#'
#'  \item{\code{Branching_times}}{Stem age of the population/species in the case of \code{"NonEndemic"},
#'  \code{"NonEndemic_MaxAge"} and \code{"Endemic"} species with no extant close relatives on
#'  the island. Set \code{NA} or \code{"NA"} if colonisation time unknown and no upper bound is known.
#'  For \code{"Endemic"} cladogenetic species these should be branching times of the
#'  radiation, including the stem age of the radiation (colonisation time estimate).
#'  If multiple times are provided in non-chronological order, they will be sorted
#'  and the oldest taken as colonisation time.}
#'
#' }
#'
#' @param island_age Age of the island in appropriate units
#'
#' @param M The size of the mainland pool, i.e the number of species that can
#' potentially colonize the island
#'
#' @param list_type2_clades If some clades are assumed to belong to a separate
#' macroevolutionary process, \code{list_type2_clades}
#' specifies the names of the clades that have a distinct macroevolutionary
#' process. The names must match those in the \code{Clade_name} column of the source
#' data table (e.g. \code{list_type2_clades = "Finches"}). If \code{NULL} (default),
#' no type 2 clades are considered.
#'
#' @param prop_type2_pool Specifies the fraction of potential mainland
#' colonists that have a distinct macroevolutionary process. Applies only if
#' there are at least some type 2 clades to consider in \code{list_type2_clades}.
#' Default to \code{NULL}, which sets the fraction to be
#' proportional to the number of clades of distinct macroevolutionary process
#' that have colonised the island. Alternatively, the user can specify a value
#' between 0 and 1 (e.g. if mainland pool size is 1000 and \code{prop_type2_pool =
#' 0.02} then number of type 2 species is 20).
#'
#' @param epss Default (\code{1e-05}) should be appropriate in most cases. This value is
#' used to set the maximum age of colonisation of \code{"NonEndemic_MaxAge"} and
#' \code{"Endemic_MaxAge"} species to an age that is slightly younger than the island
#' for cases when the age provided for that species is older than the island.
#' The new maximum age is then used as an upper bound to integrate over all
#' possible colonisation times.
#'
#' @param verbose Boolean. States if intermediate results should be printed to
#' console. Defaults to \code{TRUE}.
#'
#' @return A list object containing data.
#'
#' The first element of the list has two or three components:
#'
#' \describe{
#'
#'    \item{\code{island_age}}{the island age}
#'
#'    Then, depending on whether a distinction between species types is
#'    made, we have:
#'
#'    \item{\code{not_present}}{the number of mainland lineages that are not present on the island}
#'
#'    or:
#'
#'    \item{\code{not_present_type1}}{the number of mainland lineages of type 1 that are not present on the island}
#'    \item{\code{not_present_type2}}{the number of mainland lineages of type 2 that are not present on the island}
#'
#' }
#'
#' The following elements of the list each
#' contain information on a single colonist lineage on the island and has 5
#' components:
#'
#' \describe{
#'
#'    \item{\code{colonist_name}}{the name of the species or clade that
#'    colonized the island}
#'
#'    \item{\code{branching_times}}{island age and stem age
#'    of the population/species in the case of \code{"NonEndemic"}, \code{"NonEndemic_MaxAge"}
#'    and \code{"Endemic"} anagenetic species. For \code{"Endemic"} cladogenetic species these
#'    are island age and branching times of the radiation including the stem age
#'    of the radiation.}
#'
#'    \item{\code{stac}}{
#'
#'      the status of the colonist:
#'
#'      \itemize{
#'
#'        \item{\code{NonEndemic_MaxAge}: 1}
#'        \item{\code{Endemic}: 2}
#'        \item{\code{Endemic_and_NonEndemic}: 3}
#'        \item{\code{NonEndemic}: 4}
#'        \item{\code{Endemic_MaxAge}: 5 (if only colonisation time was given)}
#'        \item{\code{Endemic_MaxAge}: 6 (if colonisation time and cladogenesis times were given)}
#'
#'      }
#'    }
#'
#'    \item{\code{missing_species}}{number of island species that were not sampled
#'    for that particular clade (only applicable for endemic clades)}
#'
#'    \item{\code{type_1or2}}{(numeric) whether the colonist belongs to type 1 or type 2}
#'
#' }
#'
#' @author Luis M. Valente
#'
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#'
#' @export DAISIE_dataprep

# Function to prepare the data
DAISIE_dataprep = function(

    datatable,
    island_age,
    M,
    list_type2_clades = NULL,
    prop_type2_pool = NULL,
    epss = 1E-5,
    verbose = TRUE

) {

  # Overwrite some arguments if needed
  if (is.null(list_type2_clades)) {

    # That have to do with type 2 clades
    list_type2_clades <- c()
    prop_type2_pool <- NULL

  }

  # Safety checks
  check_datatable(datatable)
  check_island_age(island_age)
  check_M(M)
  check_list_type2_clades(list_type2_clades, datatable$Clade_name)
  check_prop_type2_pool(prop_type2_pool)

  # Convert to integers the numbers that must be so
  M <- as.integer(M)

  # Get the number of colonists
  n_colonisations <- nrow(datatable)

  # How many of type 2?
  n_type2_colonisations <- length(list_type2_clades)

  # The rest are type 1
  n_type1_colonisations <- n_colonisations - n_type2_colonisations

  # Type 2 mainland pool is proportional if unspecified
  if (is.null(prop_type2_pool))
    prop_type2_pool <- n_type2_colonisations / n_colonisations

  # Note: will be zero if no type 2 colonists

  # Number of clades not present
  not_present_type1 <- DDD::roundn(M * (1 - prop_type2_pool)) - n_type1_colonisations
  not_present_type2 <- DDD::roundn(M * prop_type2_pool) - n_type2_colonisations

  # Prepare empty output
  datalist <- vector("list", n_colonisations + 1L)

  # Fill in the first element
  datalist[[1]] <- list(

    island_age = island_age,
    not_present = not_present_type1

  )

  # If needed...
  if (n_type2_colonisations > 0) {

    # Add type 2 species to it
    datalist[[1]] <- c(datalist[[1]], not_present_type2 = not_present_type2)
    names(datalist[[1]])[2] <- "not_present_type1"

  }

  # For each colonist...
  for (i in 1:nrow(datatable)) {

    # Extract useful information for that clade
    clade_name <- datatable[["Clade_name"]][i]
    missing_species <- datatable[["Missing_species"]][i]
    status <- datatable[["Status"]][i]

    # Extract branching times as numbers
    brts <- get_branching_times(

      datatable[["Branching_times"]][i],
      island_age, epss, verbose, clade_name

    )

    # If colonisation time has been capped to island age...
    if (brts[2] == island_age - epss & status %in% c("Endemic", "NonEndemic")) {

      # Update endemicity status of the clade accordingly
      status <- paste0(status, "_MaxAge")

    }

    # Set the status of the colonist
    stac <- switch(

      status,
      NonEndemic_MaxAge = 1L,
      Endemic = 2L,
      Endemic_and_NonEndemic = 3L,
      NonEndemic = 4L,
      Endemic_MaxAge = ifelse(length(brts) == 2L, 5L, 6L ),
      NonEndemic_MaxAge_MinAge = 8L,
      Endemix_MaxAge_MinAge = 9L

    )

    # In case of non-endemic clade...
    if (status %in% c("NonEndemic", "NonEndemic_MaxAge")) {

      # Error if branchings other than colonisation were provided
      if (length(brts) > 2L) error_too_many_brts_for_nonendemic(clade_name)

      # Error if missing species were provided
      if (missing_species > 0L) error_too_many_missing_species_for_nonendemic(clade_name)

    }

    # Is the current clade of type 2?
    type1or2 <- ifelse(clade_name %in% list_type2_clades, 2L, 1L)

    # Record all the data needed for that clade
    datalist[[i + 1]] <- list(

      colonist_name = clade_name,
      branching_times = brts,
      stac = stac,
      missing_species = missing_species,
      type1or2 = type1or2

    )

  } # end of loop through clades

  return(datalist)

}
