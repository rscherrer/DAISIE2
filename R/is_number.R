# Function to combine with NAs
include_nas <- function(is, is_na) {

  # is: the output so far for non-NA values
  # is_na: which of the original values were NAs

  # Prepare final output
  is0 <- rep(FALSE, length(is_na))
  is0[!is_na] <- is

  return(is0)

}

# Tell if an object is a (or several) number(s) (and what kind)
is_number <- function(x, scalar = TRUE, integer = FALSE, sign = NULL, strict = FALSE) {

  # x: the object
  # scalar: whether to count only single numbers (i.e. not vectors)
  # integer: whether to say if the number boils down to an integer
  # sign: number(s) that if provided, whether x is of the same sign (e.g. -1 or c(1, -1, 0))
  # strict: whether numbers should be considered only if strictly of that sign (i.e. zero does not count as positive and negative)

  # Sanity check
  testit::assert(is_yes_no(scalar))
  testit::assert(is_yes_no(integer))
  testit::assert(is_yes_no(strict))

  # NULL or empty stuff
  if (length(x) == 0L) return(FALSE)

  # Extra checks
  if (!is.null(sign)) testit::assert(is.numeric(sign))
  if (!is.null(sign)) testit::assert(length(sign) %in% c(1L, length(x)))

  # Is the object numerical?
  if (!is.numeric(x)) return(FALSE)

  # If so, and if we want only single numbers, exit if needed
  if (scalar & length(x) > 1L) return(FALSE)

  # Find the NAs
  is_na <- is.na(x)

  # Restrict to non-NAs
  x <- x[!is_na]

  # Benefit of the doubt
  is <- rep(TRUE, length(x))

  # Only keep integer(s) if needed
  if (integer) is <- x == floor(x)

  # Exit if no need to check sign
  if (is.null(sign)) return(include_nas(is, is_na))

  # Duplicate the sign if needed
  if (length(sign) == 1L) sign <- rep(sign, length(x))

  # If non strict, count zero as both positive and negative
  if (!strict) x[x == 0] <- sign[x == 0]

  # Only keep the elements of the required sign(s)
  is <- is & sign(x) == sign

  # Exit
  return(include_nas(is, is_na))

}

# Wrappers
is_positive <- function(x, strict = FALSE) is_number(x, scalar = TRUE, sign = 1, strict)
is_negative <- function(x, strict = FALSE) is_number(x, scalar = TRUE, sign = -1, strict)
is_positive_integer <- function(x, strict = FALSE) is_number(x, scalar = TRUE, integer = TRUE, sign = 1, strict)

# Function to determine if something is a numeric vector
is_numeric_vector <- function(x, sign = NULL, strict = FALSE, na = FALSE) {

  # Is it a numeric vector?
  is <- is_number(x, scalar = FALSE, sign = sign, strict)

  # A single value is not numeric if it is NA
  if (length(is) == 1L) return(is)

  # If NAs are allowed as part of the vector...
  if (na) is[is.na(x)] <- TRUE

  return(all(is))

}

# Wrappers
is_negative_vector <- function(x, strict = FALSE, na = FALSE) is_numeric_vector(x, sign = -1, strict, na)
is_positive_vector <- function(x, strict = FALSE, na = FALSE) is_numeric_vector(x, sign = 1, strict, na)
