# Tell if an object is a (or several) number(s) (and what kind)
is_number <- function(x, scalar = TRUE, integer = FALSE, sign = NULL, strict = FALSE) {

  # x: the object
  # scalar: whether to count only single numbers (i.e. not vectors)
  # integer: whether to say if the number boils down to an integer
  # sign: number(s) that if provided, whether x is of the same sign (e.g. -1 or c(1, -1, 0))
  # strict: whether numbers should be considered only if strictly of that sign (i.e. zero does not count as positive and negative)

  # Sanity check
  testit::assert(is.logical(scalar))

  # Is the object numerical?
  if (!is.numeric(x)) return(FALSE)

  # If so, and if we want only single numbers, exit if needed
  if (scalar) if (length(x) > 1L) return(FALSE)

  # Prepare output (from here x could also be a vector)
  is <- rep(TRUE, length(x))

  # Only keep integer(s) if needed
  if (integer) is <- x == floor(x)

  # Exit if no need to check sign
  if (is.null(sign)) return(is)

  # By now sign should be a number or a vector of numbers as long as x
  testit::assert(is.numeric(sign))
  testit::assert(length(sign) %in% c(1L, length(x)))

  # If non strict, count zero as both positive and negative
  if (!strict) x[x == 0] <- sign[x == 0]

  # Only keep the elements of the required sign(s)
  is <- is & sign(x) == sign

  # Revert any NAs to being FALSE
  is[is.na(x)] <- FALSE

  # Note: those hidden in an otherwise numeric vector will have been missed.

  return(is)

}

# Wrappers
is_positive <- function(x, strict = FALSE) is_number(x, scalar = TRUE, sign = 1, strict)
is_negative <- function(x, strict = FALSE) is_number(x, scalar = TRUE, sign = -1, strict)
is_positive_integer <- function(x, strict = FALSE) is_number(x, scalar = TRUE, integer = TRUE, sign = 1, strict)
is_negative_vector <- function(x, strict = FALSE) all(is_number(x, scalar = FALSE, sign = -1, strict))
is_positive_vector <- function(x, strict = FALSE) all(is_number(x, scalar = FALSE, sign = 1, strict))
is_numeric_vector <- function(x) all(is_number(x, scalar = FALSE))
