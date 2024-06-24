# Function to find maximum likelihood estimates
daisie_ml <- function() {

  # Perform the likelihood search
  DDD::optimizer(

    # The likelihood function
    fun = integrate_daisie,

    # Vector of initial guesses for parameter values
    trparsopt = init,

    # Other arguments of the likelihood function
    island_age = island_age,
    nmax = nmax,
    is_

    idparsopt = idparsopt,
    trparsfix = trparsfix,
    idparsfix = idparsfix,
    idparsnoshift = idparsnoshift,
    idparseq = idparseq,
    pars2 = pars2,
    datalist = datalist,
    methode = methode,
    CS_version = CS_version,
    abstolint = tolint[1],
    reltolint = tolint[2],

    # Parameters of the optimizer
    optimmethod = optimmethod,
    optimpars = optimpars,
    num_cycles = num_cycles,
    jitter = jitter

  )

}
