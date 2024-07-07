# Function to create a clade object
make_clade <- function(

  is_present = FALSE, tcol = NULL, branching_times = c(),
  tmax = NULL, tmin = NULL

) {

  # is_present: whether or not the mainland relative is on the island
  # tcol: known colonization time (could be NULL)
  # branching_times: vector of branching times (could be empty)
  # tmax: maximum colonization time (could be NULL)
  # tmin: minimum colonization time (could be NULL)

  # TODO: Might bring the checks to this function.

  # If no colonization is provided

  # Assemble into a clade object
  return(list(
    is_present = is_present,
    tcol = tcol,
    branching_times = branching_times,
    tmax = tmax,
    tmin = tmin
  ))

}
