rm(list = ls())

library(DAISIE)
library(DAISIE2)

# Load the data
data("frogs_datatable")

# Prepare the data
datalist <- DAISIE_dataprep(frogs_datatable, island_age = 30, M = 300)

# Initial values of the parameters to optimize
initparsopt <- c(0.18, 0.03, 0.0006, 2)

# Indices of the parameters to optimize (here lambda_c, mu, gamma and lambda_a)
idparsopt <- c(1, 2, 4, 5)

# Type of diversity-dependence (no diversity-dependence)
ddmodel <- 0

# Values of the fixed parameters
parsfix <- Inf

# Indices of the fixed parameters (here K is the 3rd)
idparsfix <- 3

# Prepare output table
out2err <- data.frame(
  lambda_c = NA,
  mu = NA,
  K = NA,
  gamma = NA,
  lambda_a = NA,
  loglik = NA,
  df = NA,
  conv = NA
)

# ???
out2err <- invisible(out2err)

# Empty vector of indices for parameters at equilibrium (???)
idparseq <- c()

# Here we assume no equilibrium
eqmodel <- 0

# Names of model parameters
namepars <- c(
  "lambda_c",
  "mu",
  "K",
  "gamma",
  "lambda_a",
  "lambda_c2",
  "mu2",
  "K2",
  "gamma2",
  "lambda_a2",
  "prop_type2"
)

# Seems to be the indices of type 2 parameters
all_no_shift <- 6:10

# Total number of model parameters
max_idpars <- length(namepars)

# Indices of type 2 parameters not to consider different
idparsnoshift <- 6:10

# That should be all the possible indices no?
idpars <- sort(c(idparsopt, idparsfix, idparsnoshift, idparseq))

# Get the numbers of missing species for each clade
missnumspec <- unlist(lapply(datalist, function(l) { l$missing_species }))

# Rescale the initial guesses of optimised parameters (infinity becomes one)
trparsopt <- initparsopt / (1 + initparsopt)
trparsopt[which(initparsopt == Inf)] <- 1

# Same for values of fixed parameters
trparsfix <- parsfix / (1 + parsfix)
trparsfix[which(parsfix == Inf)] <- 1

# Type of conditioning
cond <- 0

# Talkative?
verbose <- 0

# No island ontogeny
island_ontogeny <- NA

# Tolerance parameters
tol <- c(1E-4, 1E-5, 1E-7)

# Maximum number of iterations
maxiter <- 1000 * round((1.25) ^ length(idparsopt))

# These two things
x_E <- 0.95
x_I <- 0.98

res <- 100

# Extra parameters
pars2 <- c(
  res,
  ddmodel,
  cond,
  verbose,
  island_ontogeny,
  eqmodel,
  tol,
  maxiter,
  x_E,
  x_I
)

# Optimisation parameters
optimpars <- c(tol, maxiter)

# Entering DAISIE:::DAISIE_loglik_all_choosepar

# Extra tolerance parameters (those are for deSolve::ode I think)
tolint <- c(1E-16, 1E-10)
abstolint <- tolint[1]
reltolint <- tolint[2]

# Non-oceanic?
non_oceanic_option <- FALSE

# Prepared another vector of transformed parameter values
trpars1 <- rep(0, 6)

# There is no type 2 clade
prop_type2_present <- which(idparsfix == 11)

# So we removed everything from the fixed parameters (REALLY?)
trparsfix <- trparsfix[-prop_type2_present]
idparsfix <- idparsfix[-prop_type2_present]

# I think there might be a problem here.

# Add values of the parameters to optimize in the vector of parameter values
trpars1[idparsopt] <- trparsopt

# Rescale again
pars1 <- trpars1 / (1 - trpars1)

# Entering DAISIE::DAISIE_loglik_all

# I think that was already the case
pars1 <- as.numeric(pars1)

# Extract this conditioning option from the extra parameters (call them options?)
cond <- pars2[3]

# That's if there are 6 values in the vector of primary parameters
endpars1 <- 6




# Entering DAISIE:::DAISIE_loglik_CS_choice to compute logp0 (wrapper)

# Entering DAISIE:::DAISIE_loglik to compute logp0

################################

# Island age
island_age <- 30

# Parameter values
lambda_c <- 0.18
mu <- 0.03
K <- 0
gamma <- 0.0006
lambda_a <- 2

# Type of diversity-dependence
ddep <- 0

# Maximum allowed number of unobserved species
max_n_unobserved <- 100

# Why is there an uneven number of state variables?



# Set up the initial probabilities
probs <- rep(0, 2 * max_n_unobserved + 1)

probs[1] <- 1 # Q^k_n
probs[max_n_unobserved + 1] <- 0 # Q^{M,k}_n

# Is this the number of observed lineages (zero for now)?
k1 <- 0

# Method for integration
methode <- "deSolve_R::lsodes"

# Entering DAISIE:::DAISIE_integrate (called on DAISIE::::DAISIE_loglik_rhs)

# Rename integrator parameters
rtol <- reltolint
atol <- abstolint
method <- methode

# Starting conditions
initprobs <- probs

# Parameters of the function to integrate
pars <- c(pars1, k1, ddep)

# Entering DAISIE::DAISIE_integrate_const (called on DAISIE::::DAISIE_loglik_rhs)

# Will probably be useless
rhs_func <- DAISIE:::DAISIE_loglik_rhs
function_as_text <- as.character(body(rhs_func)[2])
do_fun_1 <- grepl(pattern = "rhs <- 0", x = function_as_text) # this one here
do_fun_2 <- grepl(pattern = "rhs <- 1", x = function_as_text)
do_fun_3 <- grepl(pattern = "rhs <- 2", x = function_as_text)

# Number of numbers of missing species?
lx <- (length(initprobs) - 1)/2

# Entering DAISIE:::DAISIE_loglik_rhs_precomp to compute parsvec

# Unpack the parameters
lac <- pars[1]
mu <- pars[2]
K <- pars[3]
gam <- pars[4]
laa <- pars[5]
kk <- pars[6]
ddep <- pars[7]

# Something close to the number of equations?
nn <- -2:(max_n_unobserved + 2 * kk + 1)
lnn <- length(nn)
nn <- pmax(rep(0, lnn), nn)

# So here we make vectors of parameter values for each equation (???)
# There are lnn elements in those vectors, what do they correspond to?
laavec <- laa * rep(1, lnn)
lacvec <- lac * rep(1, lnn)
muvec <- mu * rep(1, lnn)
gamvec <- gam * rep(1, lnn)

# Exiting DAISIE:::DAISIE_loglik_rhs_precomp to compute parsvec

# Long vector of parameters for each equation
parsvec <- c(laavec, lacvec, muvec, gamvec, nn, kk)

# Entering DAISIE:::DAISE_ode_cs to compute y

# Corresponds to *_rhs I think
runmod <- "daisie_runmod"

# Number of state variables
N <- length(initprobs)

# Number of observed species? (zero for now)
kk <- parsvec[length(parsvec)]

# Number of numbers of missing species?
lx <- (N - 1) / 2

# RHS function
rhs_func <- DAISIE_loglik_rhs

# Update method by removing deSolve_R::
methode <- substring(methode, 12)

# Times from island age to present
times <- c(-island_age, 0)

# Solve the system of differential equations
y <- deSolve::ode(

  y = initprobs,
  times = c(-island_age, 0),
  func = rhs_func,
  parms = parsvec,
  atol = atol,
  rtol = rtol,
  method = methode

)[, 1:(N + 1)]

# Strop the first row and first column
probs <- y[-1,-1]

# Exiting DAISIE:::DAISE_ode_cs to compute y

# Integrated probabilities
y <- probs

# Exiting DAISIE::DAISIE_integrate_const

# Exiting DAISIE::DAISIE_integrate

# Retrieve the integrated probabilities
probs <- y

# Check probabilities
cp <- DAISIE:::checkprobs2(lv = 2 * lx, loglik, probs, verbose)

# Extract log-likelihood and probabilities
loglik <- cp[[1]]
probs <- cp[[2]]

# Update the log-likelihood
loglik <- loglik + log(probs[1])

# Exiting DAISIE:::DAISIE_loglik

# Exiting DAISIE:::DAISIE_loglik_CS_choice

# We should have logp0 by now
logp0 <- loglik
