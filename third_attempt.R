## Let us make a system of differential equations to integrate the process
## throughout the entire life of the island and without any lineage being
## making it to the present (k = 0).

rm(list = ls())

library(deSolve)

source("pad.R")

# Parameter values
pars <- list(lambda_c = 0.18, mu = 0.01, gamma = 0.2, lambda_a = 2)

# Time bounds (from when to when to integrate?)
from <- -30
to <- 0

# Maximum number of unobserved species allowed
nmax <- 5

# Number of numbers of unobserved species (incl. zero)
N <- nmax + 1L

# Initialize probabilities (the island starts with nothing on it, Qkn(0,0) = 1)
Qkn <- QMkn <- rep(0, N)
Qkn[1] <- 1

# Left and right padding for the present set of equations
left <- 2L
right <- 1L

# Note: this is because we sometimes need values from e.g. n-2 or n+1
# to compute the derivative at n, and this would error at the boundaries.

# Padding on both sides (for efficient treatment of boundary conditions)
Qkn <- pad(Qkn, left, right)
QMkn <- pad(QMkn, left, right)

# Note: do the padding outside of the right-hand side function because
# it is inefficient to do it many many times.

# Assemble them into a large system of state variables
states <- c(Qkn, QMkn)

# TODO: implement checks and safety nets.

# Right-hand side of the system of differential equations
right_hand_side <- function(Qkn, QMkn, lambda_c, mu, gamma, lambda_a, ii) {

  # Qkn: first set of probabilities (padded)
  # QMkn: second set of probabilities (padded)
  # lambda_c, mu, gamma, lambda_a: model parameters
  # ii: positions of true state variables (i.e. not those used for padding)

  # Initialize a (padded) vector of derivatives
  dQkn <- rep(0, length(Qkn))
  dQMkn <- rep(0, length(QMkn))

  # Vector of possible numbers of unobserved species
  n <- seq(ii) - 1L

  # Compute the first set of derivatives
  dQkn[ii] <- mu * QMkn[ii] + lambda_a * QMkn[ii - 1] + lambda_c * QMkn[ii - 2] +
    lambda_c * (n - 1) * Qkn[ii - 1] + mu * (n + 1) * Qkn[ii + 1] -
    (mu + lambda_c) * n * Qkn[ii] - gamma * Qkn[ii]

  # Compute the second set of derivatives
  dQMkn[ii] <- gamma * Qkn[ii] + lambda_c * (n - 1) * QMkn[ii - 1] +
    mu * (n + 1) * QMkn[ii + 1] - (mu + lambda_c) * n * QMkn[ii] -
    (mu + lambda_a + lambda_c) * QMkn[ii]

  # Combine
  return(c(dQkn, dQMkn))

}

# Dimensions of the system (used to identify which cell is which variable)
dim0 <- list(l = left, m = N, r = right)
dims <- list(Qkn = dim0, QMkn = dim0)

# Function to wrap the right-hand side into something acceptable by the solver
wrapped_rhs <- function(t, y, parms, dims) {

  # t: time
  # y: vector of state variables
  # parms: model parameters
  # dims: a list defining the dimensions of the system

  # Use the dimensions to identify each set of probabilities
  iQkn <- with(dims$Qkn, seq(l + m + r))
  iQMkn <- iQkn[length(iQkn)] + with(dims$QMkn, seq(l + m + r))

  # Which are the true state variables within each set?
  ii <- with(dims$Qkn, l + 1:m)

  # Pass everything to the right-hand side function
  dy <- with(pars, right_hand_side(y[iQkn], y[iQMkn], lambda_c, mu, gamma, lambda_a, ii))

  # Return the computed derivatives as a list (requirement of deSolve::ode)
  return(list(dy))

}

# Integrate
ode(y = states, times = c(from, to), func = wrapped_rhs, parms = pars, dims = dims)
