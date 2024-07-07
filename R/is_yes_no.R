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
