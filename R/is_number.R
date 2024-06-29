# Function to test whether an input is a single non-NA logical
is_yes_no <- function(x) {

  # TODO: Make use of this function wherever needed

  # Check that the input is TRUE or FALSE and nothing else
  if (!is.logical(x)) return(FALSE)
  if (length(x) > 1L) return(FALSE)
  if (is.na(x)) return(FALSE)

  # If it passes all the filters, return TRUE
  return(TRUE)

}

# Tell if an object is a (or several) number(s) (and what kind)
is_number <- function(x, scalar = TRUE, integer = FALSE, sign = NULL, strict = FALSE) {

  # x: the object
  # scalar: whether to count only single numbers (i.e. not vectors)
  # integer: whether to say if the number boils down to an integer
  # sign: number(s) that if provided, whether x is of the same sign (e.g. -1 or c(1, -1, 0))
  # strict: whether numbers should be considered only if strictly of that sign (i.e. zero does not count as positive and negative)

  # TODO: Replace sign with positive in use cases of is_number.

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

  return(is)

}
