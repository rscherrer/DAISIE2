# Constant left and right padding (same everywhere in the package)
get_left <- function() { return(2L) }
get_right <- function() { return(1L) }

# Function to get the number of variables per set (including padding)
get_nps <- function(N) {

  # N: number of possible numbers of unobserved species

  # Padd left and right
  get_left() + N + get_right()

}
