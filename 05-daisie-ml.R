# Function to find maximum likelihood estimates
daisie_ml <- function(

  data, pars, island_age, M, nmax,
  optimmethod = "subplex",
  tol = c(1e-4, 1e-5, 1e-7),
  maxiter = 1000L * round((1.25)^length(pars)),
  num_cycles = 1L, jitter = 0

) {

  # Perform the likelihood search
  DDD::optimizer(

    # The likelihood function
    fun = calc_loglik,

    # Vector of initial guesses for parameter values
    trparsopt = pars,

    # Other arguments of the likelihood function
    island_age = island_age,
    M = M,
    nmax = nmax,

    # Parameters of the optimizer
    optimmethod = optimmethod,
    optimpars = c(tol, maxiter),
    num_cycles = num_cycles,
    jitter = jitter

  )
}
