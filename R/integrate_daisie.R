# Function to return names of model parameters
par_names <- function() {

  return(c("lambda_c", "mu", "gamma", "lambda_a"))

}

# Check parameters passed to a function
check_pars <- function(pars) {

  # Check that the object is a list
  testit::assert(is.list(pars))

  # Check that all the expected parameters are there
  testit::assert(all(names(pars) %in% par_names()))
  testit::assert(all(par_names() %in% names(pars)))

  # Check that all parameters are positive numbers
  testit::assert(is_number(unlist(pars), scalar = FALSE, sign = 1))

}

# Function to check that the timings are good
check_times <- function(island_age, tcol, tmin, tmax, branching_times) {

  # island_age: age of the island
  # tcol: known colonization time (could be NULL)
  # tmin: minimum colonization time (could be NULL)
  # tmax: maximum colonization time (could be NULL)
  # branching_times: vector of branching times (could be empty)

  # Check that island age is a negative number
  testit::assert(is_number(island_age, sign = -1))

  # If known colonization time, upper and lower bounds should not be set
  if (!is.null(tcol)) testit::assert(is.null(tmax) & is.null(tmin))

  # If upper bound given, there should be no known colonization time
  if (!is.null(tmax)) testit::assert(is.null(tcol))

  # If a lower bound is given, then an upper bound should have been given too
  if (!is.null(tmin)) testit::assert(!is.null(tmax))

  # If no colonization at all, there cannot be any branching event
  if (all(unlist(lapply(list(tcol, tmin, tmax), is.null))))
    testit::assert(length(branching_times) == 0L)

  # Check that the colonization times are all negative numbers
  if (!is.null(tcol)) testit::assert(is_number(tcol, sign = -1))
  if (!is.null(tmax)) testit::assert(is_number(tmax, sign = -1))
  if (!is.null(tmin)) testit::assert(is_number(tmin, sign = -1))

  # And that if provided branching times are too
  if (length(branching_times) > 0L)
    testit::assert(all(is_number(branching_times, scalar = FALSE, sign = -1)))

  # Make sure that timings are in chronological order
  times <- c(island_age, tmax, tcol, tmin, branching_times)
  testit::assert(all(diff(times) > 0))

}

# Function to integrate the dynamics throughout island life for a given clade
integrate_daisie <- function(

  island_age, pars, nmax, tcol = NULL, branching_times = c(),
  tmax = NULL, tmin = NULL

) {

  # island_age: the age of the island
  # pars: list of model parameters
  # nmax: maximum allowed number of unobserved species
  # tcol: known colonization time
  # branching_times: vector of times of cladogenesis events
  # tmax: upper bound for an unknown colonization time
  # tmin: lower bound for an unknown colonization time

  # Check the list of parameters
  check_pars(pars)

  # Check that the maximum number of unobserved species is a positive number
  testit::assert(is_number(nmax, scalar = TRUE, integer = TRUE, sign = 1))

  # Check the provided timings
  check_times(island_age, tcol, tmin, tmax, branching_times)

  # The number of observed species we should get to at the end of the integration
  kend <- length(branching_times) + 1L * any(!is.null(c(tcol, tmin, tmax)))

  # Modify timings as needed (positive numbers so will never be reached)
  if (is.null(tmax)) tmax <- 1
  if (is.null(tmin)) tmin <- 1
  if (is.null(tcol)) tcol <- 1
  if (is.null(branching_times)) branching_times <- 1

  # Number of possible numbers of unobserved species (incl. zero)
  N <- nmax + 1L

  # Number of state variables per (padded) set of equations
  nps <- get_nps(N)

  # Vector of possible numbers of unobserved species (first of which is n = 0)
  n <- seq(N) - 1L

  # Indices of the true state variables (i.e. excl. padding) within each set
  ii <- get_left() + seq(N)

  # Index of the state variable corresponding to n = 0
  izero <- ii[1]

  # Initialize the state variables
  Qkn <- QMkn <- rep(0, nps)
  Qkn[izero] <- 1

  # Combine the state variables into a single vector
  Q <- c(Qkn, QMkn)

  # Indices of state variables for the first set
  jj <- seq(nps)

  # Identify the different sets of probabilities
  iQkn <- jj
  iQMkn <- jj + nps

  # Initialize the number of species
  k <- 0L

  # Initialize time
  t <- island_age

  # Iterator for next time point
  inext <- 2L

  # Branching time iterator
  ibrt <- 1L

  # Time points
  time_points <- c(island_age, tmax, tcol, tmin, branching_times, 0)

  # Remove those that will never be reached (those are set as positive values)
  time_points <- time_points[time_points <= 0]

  # Set absurd previous time point (will be updated, defined here to catch bugs)
  tprev <- 1

  # While we have not reached the present...
  while(t < 0) {

    # If we are at a maximum colonization time...
    if (t == tmax) {

      # The previous time step should be island age
      testit::assert(tprev == island_age)

      # We should still have zero observed species
      testit::assert(k == 0L)

      # Check that we only had two sets of equations coming in
      testit::assert(length(Q) == 2L * nps)

      # Re-initialize accordingly for the coming time interval
      QkMn <- Q[iQkn]
      QMkMn <- Q[iQMkn]
      Qkn <- rep(0, nps)
      QMkn <- rep(0, nps)

      # Re-combine the state variables together
      Q <- c(Qkn, QMkn, QkMn, QMkMn)

    }

    # If we are at a minimum colonization time...
    if (t == tmin) {

      # The previous time step should be a maximum colonization time
      testit::assert(tprev == tmax)

      # There should still be no observed species (for sure) yet
      testit::assert(k == 0L)

      # There should be four sets of equations (from max. colonization time)
      testit::assert(length(Q) == 4L * nps)

      # Re-initialize accordingly
      QkMn <- Qkn <- QMkn <- rep(0, nps)
      QkMn[izero] <- QMkn[izero]
      Qkn[ii] <- (n + 1) * Q[iQkn][ii + 1]
      QMkn[ii] <- (n + 1) * Q[iQMkn][ii + 1]

      # Re-combine the state variables together
      Q <- c(Qkn, QMkn, QkMn)

      # We now have one observed species for sure
      k <- 1L

    }

    # If we are at a known colonization time...
    if (t == tcol) {

      # The previous time point should be the island age
      testit::assert(tprev == island_age)

      # There should be no observed species before that
      testit::assert(k == 0L)

      # There should be two sets of equations before that (from island age)
      testit::assert(length(Q) == 2L * nps)

      # Re-initialize accordingly
      QkMn <- with(pars, gamma * Qkn + gamma * QMkn)
      Qkn <- Q[iQkn]
      QMkn <- Q[iQMkn]

      # Re-combine the state variables together
      Q <- c(Qkn, QMkn, QkMn)

      # We now have one observed species for sure
      k <- 1L

    }

    # If we are at a branching time...
    if (t == branching_times[ibrt]) {

      # If we come from an unknown colonization time (maximum)...
      if (tprev == tmax) {

        # Then there were zero observed species before
        testit::assert(k == 0L)

        # Re-initialize the system accordingly
        Qkn <- Q[iQkn]
        QMkn <- Q[iQMkn]
        Qkn[ii] <- with(pars, lambda_c * ((n + 1) * Qkn[ii + 1] + QMkn[ii]))
        QMkn[ii] <- with(pars, lambda_c * (n + 1) * QMkn[ii + 1])

        # And add one extra species because we go from k = 0 to k = 2
        k <- k + 1L

      } else if (tprev %in% c(tmin, tcol)) {

        # Or, if it was a minimum or known colonization...

        # Then there was one observed species before
        testit::assert(k == 1L)

        # Re-initialize differently
        Qkn <- with(pars, lambda_c * (QkMn + Qkn))
        QMkn <- with(pars, lambda_c * QMkn)

      } else {

        # Otherwise, it must have been a branching time
        testit::assert(ibrt > 1L)
        testit::assert(tprev == branching_times[ibrt - 1L])

        # Check that there were at least two observed species already
        testit::assert(k > 1L)

        # Re-initialize accordingly
        Qkn <- with(pars, lambda_c * Qkn)
        QMkn <- with(pars, lambda_c * QMkn)

      }

      # Re-combine the state variables together
      Q <- c(Qkn, QMkn)

      # We now have two observed species
      k <- k + 1L

      # Increment the branching time identifier
      ibrt <- ibrt + 1L

    }

    # Coming time point
    tnext <- time_points[inext]

    # Integrate the system from the present time to the next time point
    Q <- deSolve::ode(Q, times = c(t, tnext), func = right_hand_side, parm = pars, N = N, k = k)

    # Strip down to a final vector of probabilities
    Q <- unname(Q[-1, -1])

    # Update time
    t <- tnext

    # And move on to the next time step
    inext <- inext + 1L

    # But keep the previous time step in memory
    tprev <- time_points[inext - 2L]

  }

  # Check that we have incremented the number of observed species correctly
  testit::assert(k == kend)

  # And of course that we have reached the present
  testit::assert(t == 0)

  # Output
  return(list(Q = Q, k = k))

}
