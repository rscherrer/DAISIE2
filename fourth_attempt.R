## Let us try to integrate the dynamics from island formation (30 Mya) to the present,
## with an uncertain colonization time between 20 and 10 Mya, resulting
## in a singleton species at present.

rm(list = ls())

source("right_hand_side.R")
source("pad.R")
source("get_nps.R")

## MODEL PRAMETERS

# Parameter values
pars <- list(lambda_c = 0.18, mu = 0.01, gamma = 0.2, lambda_a = 2)


## ISLAND PARAMETERS

# Island age
island_age <- 30

# Maximum number of unobserved species allowed
nmax <- 10L

# Number of numbers of unobserved species (incl. zero)
N <- nmax + 1L

# Number of variables per set of equations (incl. padding)
nps <- get_nps(N)


## CLADE PARAMETERS

# Maximum colonization
tmax <- -20

# Minimum colonization time
tmin <- -10

# True colonization time
tc <- NULL

# Branching times
branching_times <- c(-9.99, -5, -3, -1)


# Initialize the system of differential equations
Qkn <- QMkn <- rep(0, N)
Qkn[1] <- 1

# Padding on both sides (for efficient treatment of boundary conditions)
Qkn <- pad(Qkn)
QMkn <- pad(QMkn)

# Note: this is because we sometimes need values from e.g. n-2 or n+1
# to compute the derivative at n, and this would error at the boundaries.

# Assemble them into a single system of state variables (probabilities)
Q <- c(Qkn, QMkn)


## TIME PERIOD PARAMETERS AT ISLAND AGE

# Integrate from when to when?
start <- island_age
end <- tmax

# How many observed species in that time period?
k <- 0L


# Integrate the system from start to end
Q <- deSolve::ode(Q, times = c(start, end), func = right_hand_side, parm = pars, N = N, k = k)

# Strip down to a final vector of probabilities
Q <- unname(Q[-1, -1])


## WE HAVE NOW REACHED A TMAX

start <- tmax
end <- tmin
k <- 0L


# Indices of state variables for the first set
jj <- seq(nps)

# Identify the different sets of probabilities
iQkn <- jj
iQMkn <- jj + nps

# Initialize new vectors of probabilities according to the equations
QkMn <- Q[iQkn]
QMkMn <- Q[iQMkn]
Qkn <- rep(0, length(Qkn))
QMkn <- rep(0, length(QMkn))

# Assemble the sets into one vector again
Q <- c(Qkn, QMkn, QkMn, QMkMn)

# Integrate the system from start to end
Q <- deSolve::ode(Q, times = c(start, end), func = right_hand_side, parm = pars, N = N, k = k)

# Strip down to a final vector of probabilities
Q <- unname(Q[-1, -1])


## WE HAVE NOW REACHED A TMIN (OR A TC)

start <- tmin
end <- branching_times[1]
k <- 1L


# Initialize new vectors of probabilities according to the equations
QkMn <- with(pars, gamma * Q[iQkn] + gamma * Q[iQMkn])

# Assemble the sets into one vector again
Q <- c(Q[iQkn], Q[iQMkn], QkMn)

# Integrate the system from start to end
Q <- deSolve::ode(Q, times = c(start, end), func = right_hand_side, parm = pars, N = N, k = k)

# Strip down to a final vector of probabilities
Q <- unname(Q[-1, -1])


## WE HAVE NOW REACHED THE FIRST BRANCHING TIME

start <- branching_times[1]
end <- branching_times[2]
k <- 2L


Q
