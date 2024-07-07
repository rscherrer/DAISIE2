# Function to compute the right-hand side of the system
right_hand_side <- function(t, Q, pars, N, k = 0L) {

  # t: time point (not used for the present ODEs, but still required)
  # Q: vector of state variables
  # pars: model parameters
  # N: number of numbers of unobserved species (incl. zero)
  # k: number of observed species

  # TODO: Could we pass the parameters separately to this function?

  # Number of variables per set (incl. padding)
  nps <- get_nps(N)

  # Indices of state variables for the first set
  jj <- seq(nps)

  # Extract the different sets of probabilities
  Qkn <- Q[jj]
  QMkn <- Q[jj + nps]

  # Number of sets of probabilities
  nsets <- length(Q) %/% nps

  # Extract extra sets if needed
  if (nsets > 2L) QkMn <- Q[jj + 2L * nps]
  if (nsets > 3L) QMkMn <- Q[jj + 3L * nps]

  # Initialize (padded) vectors of derivatives
  dQkn <- rep(0, nps)
  dQMkn <- rep(0, nps)
  if (nsets > 2L) dQkMn <- rep(0, nps)
  if (nsets > 3L) dQMkMn <- rep(0, nps)

  # Vector of possible numbers of unobserved species (first of which is n = 0)
  n <- seq(N) - 1L

  # Indices of the true state variables (i.e. excl. padding) within each set
  ii <- get_left() + seq(N)

  # Note: deSolve requires this function to return a list.

  # In case of a system with only two sets of equations...
  if (nsets == 2L) {

    # Compute the first set of derivatives
    dQkn[ii] <- with(pars, lambda_c * (n + 2 * k - 1) * Qkn[ii - 1] +
      lambda_a * QMkn[ii - 1] +
      lambda_c * QMkn[ii - 2] +
      mu * (n + k + 1) * Qkn[ii + 1] +
      mu * QMkn[ii] -
      (lambda_c + mu) * (n + k) * Qkn[ii] -
      gamma * Qkn[ii])

    # Compute the second set of derivatives
    dQMkn[ii] <- with(pars, gamma * Qkn[ii] +
      lambda_c * (n + 2 * k - 1) * QMkn[ii - 1] +
      mu * (n + 1) * QMkn[ii + 1] -
      (lambda_a + lambda_c + mu) * QMkn[ii] -
      (lambda_c + mu) * (n + k) * QMkn[ii])

    return(list(c(dQkn, dQMkn)))

  } else if (nsets == 3L) {

    # Otherwise, in case of three sets, compute the first set
    dQkn[ii] <- with(pars, 2 * lambda_c * QkMn[ii - 1] +
      lambda_a * QkMn[ii] +
      lambda_c * (n + 1) * Qkn[ii - 1] +
      lambda_a * QMkn[ii - 1] +
      lambda_c * QMkn[ii - 2] +
      mu * (n + 1) * Qkn[ii + 1] +
      mu * QMkn[ii] -
      (lambda_c + mu) * (n + 1) * Qkn[ii] -
      gamma * Qkn[ii])

    # Compute the second set of derivatives
    dQMkn[ii] <- with(pars, gamma * Qkn[ii] +
      lambda_c * (n + 1) * QMkn[ii - 1] +
      mu * (n + 1) * QMkn[ii + 1] -
      (lambda_a + lambda_c + mu) * QMkn[ii] -
      (lambda_c + mu) * (n + 1) * QMkn[ii])

    # Compute the third set of derivatives
    dQkMn[ii] <- with(pars, lambda_c * (n - 1) * QkMn[ii - 1] +
      mu * (n + 1) * QkMn[ii + 1] -
      (lambda_c + mu) * (n + 1) * QkMn[ii] -
      (lambda_a + gamma) * QkMn[ii])

    return(list(c(dQkn, dQMkn, dQkMn)))

  } else {

    # Or else, there are four sets: compute the first set of derivatives
    dQkn[ii] <- with(pars, lambda_c * (n - 1) * Qkn[ii - 1] +
      lambda_a * QMkn[ii - 1] +
      lambda_c * QMkn[ii - 2] +
      mu * (n + 1) * Qkn[ii + 1] +
      mu * QMkn[ii] -
      (lambda_c + mu) * n * Qkn[ii] -
      gamma * Qkn[ii])

    # Compute the second set of derivatives
    dQMkn[ii] <- with(pars, gamma * Qkn[ii] +
      gamma * QkMn[ii] +
      gamma * QMkMn[ii] +
      lambda_c * (n - 1) * QMkn[ii - 1] +
      mu * (n + 1) * QMkn[ii + 1] -
      (lambda_a + lambda_c + mu) * QMkn[ii] -
      (lambda_c + mu) * n * QMkn[ii])

    # Compute the third set of derivatives
    dQkMn[ii] <- with(pars, lambda_c * (n - 1) * QkMn[ii - 1] +
      lambda_a * QMkMn[ii - 1] +
      lambda_c * QMkMn[ii - 2] +
      mu * (n + 1) * QkMn[ii + 1] +
      mu * QMkMn[ii] -
      (lambda_c + mu) * n * QkMn[ii] -
      gamma * QkMn[ii])

    # Compute the fourth set of derivatives
    dQMkMn[ii] <- with(pars, lambda_c * (n - 1) * QMkMn[ii - 1] +
      mu * (n + 1) * QMkMn[ii + 1] -
      (lambda_a + lambda_c + mu) * QMkMn[ii] -
      (lambda_c + mu) * n * QMkMn[ii] -
      gamma * QMkMn[ii])

    return(list(c(dQkn, dQMkn, dQkMn, dQMkMn)))

  }
}
