# Function to map real numbers to a new range
range_transform <- function(x, xmax = 1, inverse = FALSE) {

  # x: vector of values to transform
  # xmax: upper bound of the transformed range (full range is between -xmax and +xmax)
  # inverse: whether to perform the inverse of the transformation

  # Check arguments
  if (!is_numeric_vector(x)) stop("x must be a numeric vector")
  if (!is_positive(xmax)) stop("xmax must be a positive number")
  if (is.infinite(xmax)) stop("xmax must be finite")
  if (!is_yes_no(inverse)) stop("inverse must be TRUE or FALSE")

  # Early exit with only zeros if needed
  if (!inverse & xmax == 0) return(rep(0, length(x)))

  # Make sure all values are within the announced range
  if (inverse & any(abs(x) > xmax)) stop("value(s) out of range detected")

  # Expand or contract the range
  y <- x * sign(x)
  y <- if (inverse) y / (xmax * sign(x) - x) else y * xmax / (sign(x) + x)

  # Edge of the original and the new space
  xlim <- ifelse(inverse, xmax, Inf)
  ylim <- ifelse(inverse, Inf, xmax)

  # Take care of edge cases
  y[x == xlim] <- ylim
  y[x == -xlim] <- -ylim
  y[x == 0] <- 0

  # Return the transformed values
  return(y)

}

# Inverse function
range_untransform <- function(x, xmax = 1) range_transform(x, xmax, inverse = TRUE)
