# Function to map real numbers to a new range
range_transform <- function(x, xmax = 1, inverse = FALSE) {

  # x: vector of values to transform
  # xmax: upper bound of the new range (full range is between -xmax and +xmax)
  # inverse: whether to perform the inverse of the transformation

  # Check arguments
  testit::assert(is_number(x, scalar = FALSE))
  testit::assert(is_number(xmax))
  testit::assert(is_yes_no(inverse))

  # If needed...
  if (inverse) {

    # Perform the inverse transformation
    y <- sign(x) * x / (xmax * sign(x) - x)

  } else {

    # Otherwise perform the regular transformation
    y <- xmax * sign(x) * x / (sign(x) + x)

  }

  # Edge of the original and the new space
  xlim <- ifelse(inverse, 1, Inf)
  ylim <- ifelse(inverse, Inf, 1)

  # Take care of edge cases
  y[x == xlim] <- ylim
  y[x == -xlim] <- -ylim
  y[x == 0] <- 0

  # Return the transformed values
  return(y)

}
