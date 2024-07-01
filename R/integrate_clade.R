# Function to guess the number of observed species
guess_k <- function(tcol, tmax, tmin, branching_times) {

  # tcol: known colonization time
  # tmax: upper bound for an unknown colonization time
  # tmin: lower bound for an unknown colonization time
  # branching_times: vector of times of cladogenesis events

  # TODO: See if we do this k business better. Is it needed?

  # If a clade has established, the size of that clade
  if (length(branching_times) > 0L) return(length(branching_times) + 1L)

  # Otherwise, if a colonization has happened for sure, k = 1
  if (!is.null(tmin) | !is.null(tcol)) return(1L)

  # Else k = 0
  return(0L)

}

# Function to pick the likelihood of the data from the set of probabilities
get_likelihood <- function(Qkn, QMkn, k, is_present, is_empty) {

  # Qkn, QMkn: sets of integrated probabilities
  # k: number of observed endemic species
  # is_present: whether the mainland species is also present
  # is_empty: did the clade leave no descendants at all?

  # Note: k does not exactly allow to distinguish all cases, as an extant
  # singleton with maximum colonization time and an empty clade both end
  # with k = 0. So use is_empty for that.

  # If the clade is absent...
  if (is_empty) return(Qkn[1])

  # Singleton where a certain (minimum) colonization time has been reached
  if (k == 1L) return(ifelse(is_present, QMkn[1], Qkn[1]))

  # Extant clade, depending on whether the mainland species is present
  if (k > 1L) return(ifelse(is_present, QMkn[1], Qkn[1]))

  # We are only left with the case of no extant species
  testit::assert(k == 0L)

  # Or because we have a singleton with maximum colonization time
  return(ifelse(is_present, QMkn[1], Qkn[2]))

}

# Function to compute the likelihood
integrate_clade <- function(

  island_age, pars, nmax, is_present = FALSE, tcol = NULL,
  branching_times = c(), tmax = NULL, tmin = NULL,
  control = list()

) {

  # island_age: the age of the island
  # pars: model parameters
  # nmax: maximum allowed number of unobserved species
  # is_present: is the mainland colonist present on the island?
  # tcol: known colonization time
  # branching_times: vector of times of cladogenesis events
  # tmax: upper bound for an unknown colonization time
  # tmin: lower bound for an unknown colonization time
  # control: named list of control parameters for deSolve::ode

  # The number of observed species we should get to at the end of the integration
  kend <- guess_k(tcol, tmax, tmin, branching_times)

  # Modify timings as needed (positive numbers so will never be reached)
  if (is.null(tmax)) tmax <- 1
  if (is.null(tmin)) tmin <- 1
  if (is.null(tcol)) tcol <- 1
  if (length(branching_times) == 0L) branching_times <- 1

  # Default control parameters if not supplied
  if (length(control) == 0L)
    control <- list(method = "lsodes", atol = 1e-16, rtol = 1e-10)

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

  # Time points
  time_points <- c(island_age, tmax, tcol, tmin, branching_times, 0)

  # Remove those that will never be reached (those are set as positive values)
  time_points <- time_points[time_points <= 0]

  # Initialize an offset we will later add to the log-likelihood
  offset <- 0

  # Note: this avoids numerical issues.

  # Initialize time
  t <- island_age

  # Iterator for next time point
  inext <- 2L

  # Branching time iterator
  ibrt <- 1L

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

      # We now have one more observed species
      k <- k + 1L

      # Increment the branching time identifier
      ibrt <- ibrt + 1L

    }

    # Coming time point
    tnext <- time_points[inext]

    # Integrate the system from the present time to the next time point
    Q <- do.call(deSolve::ode, c(

      y = list(Q),
      times = list(c(t, tnext)),
      func = right_hand_side,
      parm = list(pars),
      N = N,
      k = k,
      control

    ))

    # TODO: Avoid parameter collision here.

    # Strip down to a final vector of probabilities
    Q <- unname(Q[-1, -1])

    # Catch probabilities slightly smaller than zero
    Q[Q < 0] <- 0

    # TODO: Make sure we make the difference between numerical imprecision and
    # actual bugs.

    # Compute the sum of all probabilities
    sumQ <- sum(Q)

    # Use it to normalize the probabilities so they sum up to one
    Q <- Q / sumQ

    # TODO: Case when the sum of probabilities is zero.

    # Also use it to update the running offset for the log-likelihood
    offset <- offset + log(sumQ)

    # Note: these two last steps are to avoid numerical issues.

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

  # Unpack (and unpad) the necessary output
  Qkn <- Q[iQkn][ii]
  QMkn <- Q[iQMkn][ii]

  # Revert unspecified times to being empty (we used one as a trick)
  if (tmax == 1) tmax <- NULL
  if (tmin == 1) tmin <- NULL
  if (tcol == 1) tcol <- NULL
  if (length(branching_times) == 1L)
    if (branching_times == 1)
      branching_times <- c()

  # TODO: This is dirty.

  # Extract probability that corresponds to the likelihood
  lik <- get_likelihood(Qkn, QMkn, k, is_present, is_empty = is.null(tmax))

  # Take its logarithm
  loglik <- log(lik)

  # Add to it the offset we have prepared
  loglik <- loglik + offset

  return(loglik)

}
