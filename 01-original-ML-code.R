# replacement for DAISIE_ode_FORTRAN
# returns striped deSolve result
DAISIE_ode_cs <- function(

  initprobs,
  tvec,
  parsvec,
  atol,
  rtol,
  methode,
  runmod = "daisie_runmod"

) {

  N <- length(initprobs)
  kk <- parsvec[length(parsvec)]

  if (runmod == "daisie_runmod") {

    lx <- (N - 1) / 2
    rhs_func <- DAISIE_loglik_rhs

  } else if (runmod == "daisie_runmod1") {

    lx <- N / 4
    rhs_func <- DAISIE_loglik_rhs1

  } else if (runmod == "daisie_runmod2") {

    lx <- N / 3
    rhs_func <- DAISIE_loglik_rhs2

  }

  if (startsWith(methode, "odeint")) {

    probs <- .Call("daisie_odeint_cs", runmod, initprobs, tvec, lx, kk, parsvec[-length(parsvec)], methode, atol, rtol)

  } else if (startsWith(methode, "deSolve_R::")) {

    methode <- substring(methode, 12)

    y <- deSolve::ode(
      y = initprobs,
      times = tvec,
      func = rhs_func,
      parms = parsvec,
      atol = atol,
      rtol = rtol,
      method = methode
    )[, 1:(N + 1)]

    probs <- y[-1,-1]

  } else {

    y <- deSolve::ode(
      y = initprobs,
      parms = c(lx + 0.,kk + 0.),
      rpar = parsvec[-length(parsvec)],
      times = tvec,
      func = runmod,
      initfunc = "daisie_initmod",
      ynames = c("SV"),
      dimens = N + 2,
      nout = 1,
      outnames = c("Sum"),
      dllname = "DAISIE",
      atol = atol,
      rtol = rtol,
      method = methode
    )[,1:(N + 1)]

    probs <- y[-1,-1]  # strip 1st row and 1st column

  }

  return(probs)

}

DAISIE_loglik_rhs <- function(t, x, parsvec) {

  # x: vector of probabilities

  rhs <- 0
  kk <- parsvec[length(parsvec)]
  lx <- (length(x) - 1)/2
  lnn <- lx + 4 + 2 * kk

  laavec <- parsvec[1:lnn]
  lacvec <- parsvec[(lnn + 1):(2 * lnn)]
  muvec <- parsvec[(2 * lnn + 1):(3 * lnn)]
  gamvec <- parsvec[(3 * lnn + 1):(4 * lnn)]
  nn <- parsvec[(4 * lnn + 1):(5 * lnn)]

  xx1 = c(0,0,x[1:lx],0)
  xx2 = c(0,0,x[(lx + 1):(2 * lx)],0)
  xx3 = x[2 * lx + 1]

  nil2lx = 3:(lx + 2)

  il1 = nil2lx+kk-1
  il2 = nil2lx+kk+1
  il3 = nil2lx+kk
  il4 = nil2lx+kk-2

  in1 = nil2lx+2*kk-1
  in2 = nil2lx+1
  in3 = nil2lx+kk

  ix1 = nil2lx-1
  ix2 = nil2lx+1
  ix3 = nil2lx
  ix4 = nil2lx-2

  dx1 <- laavec[il1 + 1] * xx2[ix1] +
    lacvec[il4 + 1] * xx2[ix4] +
    muvec[il2 + 1] * xx2[ix3] +
    lacvec[il1] * nn[in1] * xx1[ix1] +
    muvec[il2] * nn[in2] * xx1[ix2] +
    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
    -gamvec[il3] * xx1[ix3]

  dx2 <- gamvec[il3] * xx1[ix3] +
    lacvec[il1 + 1] * nn[in1] * xx2[ix1] +
    muvec[il2 + 1] * nn[in2] * xx2[ix2] +
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
    -laavec[il3 + 1] * xx2[ix3]

  dx3 <- 0

  return(list(c(dx1,dx2,dx3)))

}

DAISIE_loglik_rhs1 <- function(t, x, parsvec) {

  rhs <- 1
  kk <- parsvec[length(parsvec)]
  lx <- (length(x))/4
  lnn <- lx + 4 + 2 * kk
  laavec <- parsvec[1:lnn]
  lacvec <- parsvec[(lnn + 1):(2 * lnn)]
  muvec <- parsvec[(2 * lnn + 1):(3 * lnn)]
  gamvec <- parsvec[(3 * lnn + 1):(4 * lnn)]
  nn <- parsvec[(4 * lnn + 1):(5 * lnn)]

  xx1 <- c(0,0,x[1:lx],0)
  xx2 <- c(0,0,x[(lx + 1):(2 * lx)],0)
  xx3 <- c(0,0,x[(2 * lx + 1):(3 * lx)],0)
  xx4 <- c(0,0,x[(3 * lx + 1):(4 * lx)],0)

  nil2lx <- 3:(lx + 2)

  il1 <- nil2lx+kk-1
  il2 <- nil2lx+kk+1
  il3 <- nil2lx+kk
  il4 <- nil2lx+kk-2

  in1 <- nil2lx+2*kk-1
  in2 <- nil2lx+1
  in3 <- nil2lx+kk
  in4 <- nil2lx-1

  ix1 <- nil2lx-1
  ix2 <- nil2lx+1
  ix3 <- nil2lx
  ix4 <- nil2lx-2

  dx1 <- lacvec[il1] * nn[in1] * xx1[ix1] +
    laavec[il1 + 1] * xx2[ix1] +
    lacvec[il4 + 1] * xx2[ix4] +
    muvec[il2] * nn[in2] * xx1[ix2] +
    muvec[il3 + 1] * xx2[ix3] +
    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
    -gamvec[il3] * xx1[ix3]

  dx2 <- gamvec[il3] * xx1[ix3] +
    gamvec[il3] * xx3[ix3] +
    gamvec[il3 + 1] * xx4[ix3] +
    lacvec[il1 + 1] * nn[in1] * xx2[ix1] +
    muvec[il2 + 1] * nn[in2] * xx2[ix2] +
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
    -laavec[il3 + 1] * xx2[ix3]

  dx3 <- lacvec[il1] * nn[in1] * xx3[ix1] +
    laavec[il1 + 1] * xx4[ix1] +
    lacvec[il4 + 1] * xx4[ix4] +
    muvec[il2] * nn[in2] * xx3[ix2] +
    muvec[il3 + 1] * xx4[ix3] +
    -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
    -gamvec[il3] * xx3[ix3]

  dx4 <- lacvec[il1 + 1] * nn[in1] * xx4[ix1] +
    muvec[il2 + 1] * nn[in2] * xx4[ix2] +
    -(lacvec[il3 + 1] + muvec[il3 + 1]) * nn[in3 + 1] * xx4[ix3] +
    -gamvec[il3 + 1] * xx4[ix3]

  return(list(c(dx1,dx2,dx3,dx4)))

}

DAISIE_loglik_rhs2 <- function(t, x, parsvec) {

  rhs <- 2
  kk <- parsvec[length(parsvec)]
  lx <- (length(x))/3
  lnn <- lx + 4 + 2 * kk
  laavec <- parsvec[1:lnn]
  lacvec <- parsvec[(lnn + 1):(2 * lnn)]
  muvec <- parsvec[(2 * lnn + 1):(3 * lnn)]
  gamvec <- parsvec[(3 * lnn + 1):(4 * lnn)]
  nn <- parsvec[(4 * lnn + 1):(5 * lnn)]

  xx1 = c(0,0,x[1:lx],0)
  xx2 = c(0,0,x[(lx + 1):(2 * lx)],0)
  xx3 = c(0,0,x[(2 * lx + 1):(3 * lx)],0)

  nil2lx = 3:(lx + 2)

  il1 = nil2lx+kk-1
  il2 = nil2lx+kk+1
  il3 = nil2lx+kk
  il4 = nil2lx+kk-2

  in1 = nil2lx+2*kk-1
  in2 = nil2lx+1
  in3 = nil2lx+kk
  in4 = nil2lx-1

  ix1 = nil2lx-1
  ix2 = nil2lx+1
  ix3 = nil2lx
  ix4 = nil2lx-2

  # inflow:
  # anagenesis in colonist when k = 1: Q_M,n -> Q^1_n; n+k species present
  # cladogenesis in colonist when k = 1: Q_M,n-1 -> Q^1_n;
  # n+k-1 species present; rate twice
  # anagenesis of reimmigrant: Q^M,k_n-1 -> Q^k,n; n+k-1+1 species present
  # cladogenesis of reimmigrant: Q^M,k_n-2 -> Q^k,n;
  # n+k-2+1 species present; rate once
  # extinction of reimmigrant: Q^M,k_n -> Q^k,n; n+k+1 species present
  # cladogenesis in one of the n+k-1 species: Q^k_n-1 -> Q^k_n;
  # n+k-1 species present; rate twice for k species
  # extinction in one of the n+1 species: Q^k_n+1 -> Q^k_n; n+k+1 species
  # present
  # outflow:
  # all events with n+k species present

  dx1 = (laavec[il3] * xx3[ix3] + 2 * lacvec[il1] * xx3[ix1]) * (kk == 1) +
    laavec[il1 + 1] * xx2[ix1] +
    lacvec[il4 + 1] * xx2[ix4] +
    muvec[il2 + 1] * xx2[ix3] +
    lacvec[il1] * nn[in1] * xx1[ix1] +
    muvec[il2] * nn[in2] * xx1[ix2] +
    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
    -gamvec[il3] * xx1[ix3]

  # inflow:
  # immigration when there are n+k species: Q^k,n -> Q^M,k_n;
  # n+k species present
  # cladogenesis in n+k-1 species: Q^M,k_n-1 -> Q^M,k_n;
  # n+k-1+1 species present; rate twice for k species
  # extinction in n+1 species: Q^M,k_n+1 -> Q^M,k_n; n+k+1+1 species present
  # outflow:
  # all events with n+k+1 species present

  dx2 <- gamvec[il3] * xx1[ix3] +
    lacvec[il1 + 1] * nn[in1] * xx2[ix1] +
    muvec[il2 + 1] * nn[in2] * xx2[ix2] +
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
    -laavec[il3 + 1] * xx2[ix3]

  # only when k = 1
  # inflow:
  # cladogenesis in one of the n-1 species: Q_M,n-1 -> Q_M,n;
  # n+k-1 species present; rate once
  # extinction in one of the n+1 species: Q_M,n+1 -> Q_M,n;
  # n+k+1 species present
  # outflow:
  # all events with n+k species present

  dx3 <- lacvec[il1] * nn[in4] * xx3[ix1] +
    muvec[il2] * nn[in2] * xx3[ix2] +
    -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
    -(laavec[il3] + gamvec[il3]) * xx3[ix3]

  return(list(c(dx1,dx2,dx3)))

}

DAISIE_loglik_rhs_precomp <- function(pars, lx) {

  # lx might be the total number of equations
  # kk might be the number of observed lineages

  # Unpack the parameters
  lac = pars[1]
  mu = pars[2]
  K = pars[3]
  gam = pars[4]
  laa = pars[5]
  kk = pars[6]
  ddep = pars[7]

  # No idea what that might be
  nn = -2:(lx + 2 * kk + 1)
  lnn = length(nn)
  nn = pmax(rep(0, lnn), nn)

  # Stuff changes based on the type of diversity-dependence

  if (ddep == 0) {

    # So here we make vectors of parameter values for each equation (???)
    # There are lnn elements in those vectors, what do they correspond to?

    laavec <- laa * rep(1, lnn)
    lacvec <- lac * rep(1, lnn)
    muvec <- mu * rep(1, lnn)
    gamvec <- gam * rep(1, lnn)

  } else if (ddep == 1) {

    laavec <- laa * rep(1, lnn)
    lacvec <- pmax(rep(0, lnn), lac * (1 - nn/K))
    muvec <- mu * rep(1, lnn)
    gamvec <- gam * rep(1, lnn)

  } else if (ddep == 2) {

    laavec <- laa * rep(1, lnn)
    lacvec <- pmax(rep(0,lnn), lac * exp(-nn/K))
    muvec <- mu * rep(1, lnn)
    gamvec <- gam * rep(1, lnn)

  } else if (ddep == 11) {

    laavec = laa * rep(1, lnn)

    #lacvec = pmax(rep(0,lnn),lac * (1 - nn/K))

    lacvec <- get_clado_rate_per_capita(
      lac = lac,
      d = 0,
      num_spec = nn,
      K = K,
      A = 1
    )

    muvec = mu * rep(1, lnn)
    muvec <- rep(1, lnn) * get_ext_rate_per_capita(mu = mu, x = 0)

    #gamvec = pmax(rep(0,lnn),gam * (1 - nn/K))

    gamvec <- get_immig_rate_per_capita(
      gam = gam,
      num_spec = nn,
      K = K,
      A = 1
    )

  } else if(ddep == 21) {

    laavec <- laa * rep(1, lnn)
    lacvec <- pmax(rep(0, lnn), lac * exp(-nn/K))
    muvec <- mu * rep(1, lnn)
    gamvec <- pmax(rep(0, lnn), gam * exp(-nn/K))

  } else if(ddep == 3) {

    laavec <- laa * rep(1, lnn)
    lacvec <- lac * rep(1, lnn)
    muvec <- mu * (1 + nn/K)
    gamvec <- gam * rep(1, lnn)

  }

  # cat("lacvec ", lacvec, "\n")
  # cat("muvec ", muvec, "\n")
  # cat("gamvec ", gamvec, "\n")
  # cat("laavec ", laavec, "\n")
  # cat("nn ", nn, "\n")
  # cat("kk ", kk, "\n")
  # cat("lx ", lx, "\n")
  # cat("lnn ", lnn)

  return(c(laavec, lacvec, muvec, gamvec, nn, kk))

}

DAISIE_integrate_const <- function(

  initprobs,
  tvec,
  rhs_func,
  pars,
  rtol,
  atol,
  method

) {

  # During code coverage, 'function_as_text' may become:
  #
  # if (TRUE) {
  #   covr::count("DAISIE_loglik_CS.R:58:3:58:25:3:25:4657:4657")
  #   lx <- (length(x) - 1)/2
  # }
  #
  # It is the 'lx <- [something]' part that we are interested in
  #
  # Use a regular expression to extract if the part that we are interested
  # in is present

  function_as_text <- as.character(body(rhs_func)[2])

  do_fun_1 <- grepl(pattern = "rhs <- 0", x = function_as_text)
  do_fun_2 <- grepl(pattern = "rhs <- 1", x = function_as_text)
  do_fun_3 <- grepl(pattern = "rhs <- 2", x = function_as_text)

  if (do_fun_1) {

    lx <- (length(initprobs) - 1)/2
    parsvec <- c(DAISIE_loglik_rhs_precomp(pars,lx))

    y <- DAISIE_ode_cs(
      initprobs,
      tvec,
      parsvec,
      atol,
      rtol,
      method,
      runmod = "daisie_runmod"
    )

    #y <- deSolve::ode(
    #    y = initprobs,
    #    times = tvec,
    #    func = DAISIE_loglik_rhs1,
    #    parms = parsvec,
    #    rtol = rtol,
    #    atol = atol,
    #    method = method
    #  )[2, -1]

  } else if (do_fun_2) {

    lx <- (length(initprobs))/4
    parsvec <- c(DAISIE_loglik_rhs_precomp(pars,lx))

    y <- DAISIE_ode_cs(
      initprobs,
      tvec,
      parsvec,
      atol,
      rtol,
      method,
      runmod = "daisie_runmod1"
    )

  } else if (do_fun_3) {

    lx <- (length(initprobs))/3
    parsvec <- c(DAISIE_loglik_rhs_precomp(pars,lx))

    y <- DAISIE_ode_cs(
      initprobs,
      tvec,
      parsvec,
      atol,
      rtol,
      method,
      runmod = "daisie_runmod2"
    )

    #y <- deSolve::ode(
    #  y = initprobs,
    #  times = tvec,
    #  func = DAISIE_loglik_rhs2,
    #  parms = parsvec,
    #  rtol = rtol,
    #  atol = atol,
    #  method = method
    #)[2, -1]

  } else {

    stop(
      "The integrand function is written incorrectly. ",
      "Value of 'function_as_text':", function_as_text
    )

  }

  return(y)

}

DAISIE_integrate <- function(

  initprobs,
  tvec,
  rhs_func,
  pars,
  rtol,
  atol,
  method

) {

  if (length(pars) <= 7) {

    # TODO COMPARE THESE 2. RHS1 SEEMS OK

    return(
      DAISIE_integrate_const(
        initprobs,
        tvec,
        rhs_func,
        pars,
        rtol,
        atol,
        method
      )
    )

  } else {

    return(
      DAISIE_integrate_time(
        initprobs,
        tvec,
        rhs_func,
        pars,
        rtol,
        atol,
        method
      )
    )

  }
}

DAISIE_loglik <- function(

  pars1,
  pars2,
  brts,
  stac,
  missnumspec,
  methode = "lsodes",
  abstolint = 1E-16,
  reltolint = 1E-10,
  verbose

) {

  # stac = status of the clade formed by the immigrant

  #  . stac == 1 : immigrant is present but has not formed an extant clade
  #  . stac == 2 : immigrant is not present but has formed an extant clade
  #  . stac == 3 : immigrant is present and has formed an extant clade
  #  . stac == 4 : immigrant is present but has not formed an extant clade,
  #  and it is known when it immigrated.
  #  . stac == 5 : immigrant is not present and has not formed an extant clade,
  #  but only an endemic species

  #  . stac == 6 : like 2, but with max colonization time
  #  . stac == 7 : like 3, but with max colonization time
  #  . stac == 8 : like 1, but with min colonization time
  #  . stac == 9 : like 5, but with min colonization time

  # warn if laa becomes Inf
  if (any(is.infinite(pars1))) {

    if (verbose) {

      message('One of the parameters is infinite.')

    }
  }

  if (is.na(pars2[4])) {

    pars2[4] = 0

  }

  ddep <- pars2[2]
  K <- pars1[3]

  # If ontogeny, lx should be bigger, i.e., modified by max area
  # if (!is.na(pars2[5])) {
  #   K <- K * pars1[9]
  # }

  if(length(pars1) == 6) {

    probability_of_init_presence <- pars1[6]
    pars1 <- pars1[-6]

  } else {

    probability_of_init_presence <- 0

  }

  # Ah, so brts can be negative, they will be sorted
  brts <- -sort(abs(as.numeric(brts)), decreasing = TRUE)

  # Problem
  if(length(brts) == 1 & sum(brts == 0) == 1) {

    stop('The branching times contain only a 0. This means the island emerged at the present which is not allowed.');
    loglik = -Inf
    return(loglik)

  }

  if (sum(brts == 0) == 0) {

    brts[length(brts) + 1] <- 0

  }

  # for stac = 0, brts will contain origin of island and 0; length = 2;
  # no. species should be 0
  # for stac = 1, brts will contain origin of island, maximum colonization time
  # (usually island age) and 0; length = 3; no. species should be 1
  # for stac = 2, brts will contain origin of island, colonization event,
  # branching times, 0; no. species should be no. branching times + 1
  # for stac = 3, brts will contain origin of island, colonization event,
  # branching times, 0; no. species should be no. branching times + 2
  # for stac = 4, brts will contain origin of island, colonization event and 0;
  # length = 3; no. species should be 1
  # for stac = 5, brts will contain origin of island, maximum colonization time
  # (usually island age), and 0; length = 2; number of species should be 1 (+ missing species)
  # for stac = 6, brts will contain origin of island, maximum colonization time
  # (usually island age), branching times and 0;
  # number of species should be no. branching times + 1
  # for stac = 7, brts will contain origin of island, maximum colonization time
  #  usually island age), branching times and 0;
  #  number of species should be no. branching times + 2
  # for stac = 8, brts will contain origin of island, maximum colonization time
  #  usually island age), minimum colonization time and 0; length = 4;
  #  number of species should be 1
  # for stac = 9, brts will contain origin of island, maximum colonization time
  #  usually island age), minimum colonization time and 0; length = 4;
  #  number of species should be 1



  # What's that thing?

  S <- 0 * (stac == 0) +
    (stac == 1 || stac == 4 || stac == 5 || stac == 8 || stac == 9) +
    (length(brts) - 2) * (stac == 2) + (length(brts) - 1) * (stac == 3) +
    (length(brts) - 2) * (stac == 6) + (length(brts) - 1) * (stac == 7)

  S2 <- S - (stac == 1) - (stac == 3) - (stac == 4) - (stac == 7)

  loglik <- -lgamma(S2 + missnumspec + 1) +
    lgamma(S2 + 1) + lgamma(missnumspec + 1)




  if (min(pars1) < 0) {

    message("One or more parameters are negative.")
    loglik <- -Inf
    return(loglik)

  }

  if ((ddep == 1 | ddep == 11) & ceiling(K) < (S + missnumspec)) {

    if (verbose) {

      message(
        'The proposed value of K is incompatible with the number of species
          in the clade. Likelihood for this parameter set
          will be set to -Inf. \n'
      )
    }

    loglik <- -Inf
    return(loglik)

  }

  lac <- pars1[1]

  if (lac == Inf & missnumspec == 0 & length(pars1) == 5) {

    if (verbose) warning('Infinite lambda detected')
    loglik <- DAISIE::DAISIE_loglik_high_lambda(pars1, -brts, stac)

  } else {

    if (ddep == 1 | ddep == 11) {

      lx <- min(
        1 + max(missnumspec, ceiling(K)),
      )

    } else {

      lx <- DDD::roundn(pars2[1]) + missnumspec

    }

    if (loglik > -Inf) {

      # in all cases we integrate from the origin of the island to the colonization event
      # (stac 2, 3, 4), the first branching point (stac = 6, 7), to the maximum colonization
      # time (stac = 1, 5, 8, 9) or to the present (stac = 0)

      probs <- rep(0,2 * lx + 1)
      probs[1] <- 1 - probability_of_init_presence #Q^k_n
      probs[lx + 1] <- probability_of_init_presence #Q^{M,k}_n

      k1 <- 0 # pars2 (containing island_ontogeny etc, should go into pars1?)

      probs <- DAISIE:::DAISIE_integrate(
        probs,
        brts[1:2],
        DAISIE_loglik_rhs,
        c(pars1, k1, ddep),
        rtol = reltolint,
        atol = abstolint,
        method = methode
      )

      cp <- DAISIE::checkprobs2(lv = 2 * lx, loglik, probs, verbose)
      loglik <- cp[[1]]
      probs <- cp[[2]]

      if (stac == 0) {

        # for stac = 0, the integration was from the origin of the island until
        # the present so we can immediately evaluate the probability of no clade
        # being present and no immigrant species.

        loglik <- loglik + log(probs[1])

      } else {

        if (stac %in% c(1, 5:9)) {

          # for stac = 1, we now integrate from the maximum colonization time
          # (usually the island age + tiny time unit) until the present, where
          # we set all probabilities where the immigrant is already present to 0
          # and we evaluate the probability of the immigrant species being
          # present, but there can be missing species.
          # for stac = 5, we do exactly the same, but we evaluate the
          # probability of an endemic species being present alone.
          # for stac = 6 and 7, integration is from the maximum colonization
          # time until the first branching time. This is the same as we did for
          # stac = 1, 5.
          # for stac = 8 and 9, integration is from the maximum colonization
          # time until the minimum colonization time.
          # In all cases we are dealing with a maximum colonization time which
          # means that any colonization that took place before this maximum
          # colonization time (including presence in the non-oceanic scenario)
          # does not count and should be followed by another colonization.
          # To allow this we introduce a third set of equations for the
          # probability that colonization might have happened before but
          # recolonization has not taken place yet (Q_M,n).

          epss <- 1.01E-5 # We're taking the risk

          if (abs(brts[2] - brts[1]) >= epss) {

            probs[(2 * lx + 1):(4 * lx)] <- probs[1:(2 * lx)]
            probs[1:(2 * lx)] <- 0

          } else {

            probs[(2 * lx + 1):(4 * lx)] <- 0

          }

          probs <- DAISIE_integrate(
            probs,
            brts[2:3],
            DAISIE_loglik_rhs1,
            c(pars1, k1, ddep),
            rtol = reltolint,
            atol = abstolint,
            method = methode
          )

          cp <- DAISIE::checkprobs2(lx, loglik, probs, verbose)
          loglik <- cp[[1]]
          probs <- cp[[2]]

          if (stac %in% c(1, 5)) {

            loglik <- loglik + log(probs[(stac == 1) * lx + (stac == 5) + 1 + missnumspec])

          } else if (stac %in% c(6, 7, 8, 9)) {

            probs2 <- rep(0, 3 * lx)
            probs2[1:(lx - 1)] <- (1:(lx - 1)) * probs[2:lx]
            probs2[(lx + 1):(2 * lx - 1)] <- (1:(lx - 1)) * probs[(lx + 2):(2 * lx)]
            probs2[2 * lx + 1] <- probs[(lx + 1)]
            probs2[(2 * lx + 2):(3 * lx)] <- 0

            probs <- probs2
            rm(probs2)

            if (stac %in% c(8, 9)) {

              k1 <- 1

              probs <- DAISIE_integrate(
                probs,
                c(brts[3:4]),
                DAISIE_loglik_rhs2,
                c(pars1, k1, ddep),
                rtol = reltolint,
                atol = abstolint,
                method = methode
              )

              cp <- checkprobs2(lx, loglik, probs, verbose)
              loglik <- cp[[1]]
              probs <- cp[[2]]

              loglik <- loglik + log(probs[(stac == 8) * (2 * lx + 1) + (stac == 9) + missnumspec])

            }
          }

        } else if (stac %in% c(2, 3, 4)) {

          # for stac = 2, 3, 4, integration is then from the colonization
          # event until the first branching time (stac = 2 and 3) or the present
          # (stac = 4). We add a set of equations for Q_M,n, the probability
          # that the process is compatible with the data, and speciation has not
          # happened; during this time immigration is not allowed because it
          # would alter the colonization time.

          t <- brts[2]

          gamvec <- divdepvec(
            lac_or_gam = "gam",
            pars1 = pars1,
            t = t,
            lx = lx,
            k1 = k1,
            ddep = ddep * (ddep == 11 | ddep == 21)
          )

          probs[(2 * lx + 1):(3 * lx)] <- gamvec[1:lx] * probs[1:lx] +
            gamvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]

          probs[1:(2 * lx)] <- 0

          k1 <- 1

          probs <- DAISIE_integrate(
            probs,
            c(brts[2:3]),
            DAISIE_loglik_rhs2,
            c(pars1, k1, ddep),
            rtol = reltolint,
            atol = abstolint,
            method = methode
          )

          cp <- checkprobs2(lx, loglik, probs, verbose)
          loglik <- cp[[1]]
          probs <- cp[[2]]

          # if stac = 4, we're done and we take an element from Q_M,n
          if (stac == 4) {

            loglik = loglik + log(probs[2 * lx + 1 + missnumspec])

          }

        } else if (stac %in% c(2, 3, 6, 7)) {

          # at the first branching point all probabilities of states Q_M,n are
          # transferred to probabilities where only endemics are present. Then
          # go through the branching points.

          S1 <- length(brts) - 1
          startk <- 3

          if (S1 >= startk) {

            t <- brts[startk]

            lacvec <- divdepvec(
              lac_or_gam = "lac",
              pars1 = pars1,
              t = t,
              lx = lx + stac %in% c(6,7),
              k1 = k1,
              ddep = ddep
            )

            if (stac %in% c(2,3)) {

              probs[1:lx] <- lacvec[1:lx] * (probs[1:lx] + probs[(2 * lx + 1):(3 * lx)])
              probs[(lx + 1):(2 * lx)] <- lacvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]

            } else {

              # stac in c(6,7)
              probs[1:(lx - 1)] <- lacvec[2:lx] *
                ((1:(lx - 1)) * probs[2:lx] + probs[(2 * lx + 1):(3 * lx - 1)])
              probs[(lx + 1):(2 * lx - 1)] <- lacvec[3:(lx + 1)] * (1:(lx - 1)) *
                probs[(lx + 2):(2 * lx)]
              probs[lx] <- 0
              probs[2 * lx] <- 0

            }

            probs <- probs[-c((2 * lx + 2):(3 * lx))]
            probs[2 * lx + 1] <- 0

            # After speciation, colonization is allowed again (re-immigration)
            # all probabilities of states with the immigrant present are set to
            # zero and all probabilities of states with endemics present are
            # transported to the state with the colonist present waiting for
            # speciation to happen. We also multiply by the (possibly diversity-
            # dependent) immigration rate.

            for (k in startk:S1) {

              k1 <- k - 1

              probs <- DAISIE_integrate(
                probs,
                brts[k:(k+1)],
                DAISIE_loglik_rhs,
                c(pars1, k1, ddep),
                rtol = reltolint,
                atol = abstolint,
                method = methode
              )

              cp <- DAISIE::checkprobs2(lx, loglik, probs, verbose)
              loglik = cp[[1]]
              probs = cp[[2]]

              if (k < S1) {

                # speciation event
                t <- brts[k + 1]

                lacvec <- divdepvec(
                  lac_or_gam = "lac",
                  pars1 = pars1,
                  t = t,
                  lx = lx,
                  k1 = k1,
                  ddep = ddep
                )

                probs[1:(2 * lx)] <- c(lacvec[1:lx], lacvec[2:(lx + 1)]) *
                  probs[1:(2 * lx)]

              }
            }
          }

          # we evaluate the probability of the phylogeny with any missing species at the present without (stac = 2 or stac = 6) or with (stac = 3 or stac = 7) the immigrant species
          loglik <- loglik + log(probs[(stac %in% c(3, 7)) * lx + 1 + missnumspec])

        }
      }
    }
  }

  if (length(pars1) == 11) {

    # CHANGE

    DAISIE::print_parameters_and_loglik(
      pars = c(stac, pars1[5:10]), # should this be 6:10, or 6:11?
      loglik = loglik,
      verbose = pars2[4],
      type = 'clade_loglik'
    )

  } else {

    print_parameters_and_loglik(
      pars = c(stac, pars1[1:5]),
      loglik = loglik,
      verbose = pars2[4],
      type = 'clade_loglik'
    )
  }

  if (is.na(loglik)) {

    message("NA in loglik encountered. Changing to -Inf.")
    loglik <- -Inf

  }

  loglik <- as.numeric(loglik)

  #testit::assert(is.numeric(loglik))

  return(loglik)

}

DAISIE_loglik_CS_choice <- function(

  pars1,
  pars2,
  datalist = NULL,
  brts,
  stac,
  missnumspec,
  methode = "lsodes",
  CS_version = 1,
  abstolint = 1E-16,
  reltolint = 1E-10,
  verbose = FALSE

) {

  if (CS_version[[1]] == 1) {

    loglik <- DAISIE_loglik(
      pars1 = pars1,
      pars2 = pars2,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose
    )

  } else if (CS_version[[1]] == 2) {

    loglik <- DAISIE_loglik_integrate(
      pars1 = pars1,
      pars2 = pars2,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose
    )

  } else if (CS_version[[1]] == 0) {

    loglik <- DAISIE_loglik_IW_M1(
      pars1 = pars1,
      pars2 = pars2,
      datalist = datalist,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose
    )

  }

  return(loglik)

}

DAISIE_loglik_all <- function(

  pars1,
  pars2,
  datalist,
  methode = "lsodes",
  CS_version = 1,
  abstolint = 1E-16,
  reltolint = 1E-10

) {

  # If island ontogeny there should be 14 parameters
  if (length(pars1) == 14) {

    if (datalist[[1]]$island_age > pars1[11]) {

      stop(
        "The island age in the area parameters is inconsistent with the island
        data."
      )
    }

    peak <- calc_peak(

      total_time = datalist[[1]]$island_age,
      area_pars = create_area_pars(
        max_area = pars1[8],
        current_area = pars1[9],
        proportional_peak_t = pars1[10],
        total_island_age = pars1[11],
        sea_level_amplitude = pars1[12],
        sea_level_frequency = pars1[13],
        island_gradient_angle = pars1[14]
      )

    )

    pars1 <- c(
      pars1,
      island_ontogeny = pars2[5],
      sea_level = pars2[6],
      datalist[[1]]$island_age,
      peak
    )
  }

  # I think that was already the case
  pars1 <- as.numeric(pars1)

  # Extract this conditioning option from the extra parameters (call them options?)
  cond <- pars2[3]

  # Ahem
  if (length(pars1) == 6) {

    endpars1 <- 6

  } else {

    endpars1 <- 5

  }

  # pars2[5] is ontogeny
  if(length(pars1) %in% c(5,6) | !is.na(pars2[5])) {

    # If ontogeny
    if(!is.na(pars2[5])) {

      endpars1 <- length(pars1)

    }

    logp0 <- DAISIE:::DAISIE_loglik_CS_choice(

      pars1 = pars1,
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint

    )

    # Some special case where we have to approximate
    if(logp0 >= 0 & pars1[2]/pars1[1] > 100) {

      logp0 <- DAISIE::approximate_logp0(
        gamma = pars1[4],
        mu = pars1[2],
        t = datalist[[1]]$island_age
      )
    }

    # Niet goed
    if(logp0 >= 0) {

      message('Positive values of loglik encountered without possibility for approximation. Setting loglik to -Inf.')

      loglik <- -Inf

      print_parameters_and_loglik(
        pars = pars,
        loglik = loglik,
        verbose = pars2[4],
        type = 'island_loglik'
      )

      return(loglik)

    }

    # As many logp0 as there are missing species

    if (is.null(datalist[[1]]$not_present)) {

      # Type 2 stuff
      loglik <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0
      numimm <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1

    } else {

      # Typ 1 stuff only
      loglik <- datalist[[1]]$not_present * logp0

      # Tha's the mainand pool
      numimm <- datalist[[1]]$not_present + length(datalist) - 1

    }

    logcond <- DAISIE:::logcondprob(numcolmin = cond, numimm = numimm, logp0 = logp0)

    # If there is at east one colonist...
    if (length(datalist) > 1) {

      # For each colonist...
      for (i in 2:length(datalist)) {

        # Revert to type 1
        datalist[[i]]$type1or2 <- 1

      }
    }

  } else {

    # Not here

    numimm <- datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2 + length(datalist) - 1
    numimm_type2 <- length(which(unlist(datalist)[which(names(unlist(datalist)) == "type1or2")] == 2))
    numimm_type1 <- length(datalist) - 1 - numimm_type2

    if (is.na(pars1[11]) == FALSE && length(pars1) == 11) {

      if (pars1[11] < numimm_type2 / numimm | pars1[11] > (1 - numimm_type1 / numimm)) {

        return(-Inf)

      }

      datalist[[1]]$not_present_type2 <- max(0, round(pars1[11] * numimm) - numimm_type2)
      datalist[[1]]$not_present_type1 <- numimm - (length(datalist) - 1) - datalist[[1]]$not_present_type2

    }

    logp0_type1 <- DAISIE_loglik_CS_choice(

      pars1 = pars1[1:5],
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint

    )

    if(logp0_type1 >= 0 & pars1[2]/pars1[1] > 100) {

      logp0_type1 <- approximate_logp0(gamma = pars1[4], mu = pars1[2], t = datalist[[1]]$island_age)

    }

    logp0_type2 <- DAISIE_loglik_CS_choice(

      pars1 = pars1[6:10],
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint

    )

    if (logp0_type2 >= 0 & pars1[7] / pars1[6] > 100) {

      logp0_type2 <- approximate_logp0(gamma = pars1[9], mu = pars1[7], t = datalist[[1]]$island_age)

    }

    loglik <- datalist[[1]]$not_present_type1 * logp0_type1 + datalist[[1]]$not_present_type2 * logp0_type2

    #logcond <- (cond == 1) *
    #  log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) *
    #                logp0_type1 +
    #                (datalist[[1]]$not_present_type2 + numimm_type2) *
    #                logp0_type2))

    logcond <- logcondprob(
      numcolmin = cond,
      numimm = c(datalist[[1]]$not_present_type1 + numimm_type1,datalist[[1]]$not_present_type2 + numimm_type2),
      logp0 = c(logp0_type1,logp0_type2)
    )

  }

  loglik <- loglik - logcond

  # Again, if multiple colonists...
  if (length(datalist) > 1) {

    #  For each one...
    for (i in 2:length(datalist)) {

      # If type 1
      if (datalist[[i]]$type1or2 == 1) {

        # Restrict parameters to type 1 parameters
        pars <- pars1[1:endpars1]

      } else {

        # Or type 2 otherwise
        pars <- pars1[6:10]

      }

      # Compute likelihood for that clade and add it to what we already have
      loglik <- loglik +
        DAISIE:::DAISIE_loglik_CS_choice(
          pars1 = pars,
          pars2 = pars2,
          datalist = datalist[[i]],
          brts = datalist[[i]]$branching_times,
          stac = datalist[[i]]$stac,
          missnumspec = datalist[[i]]$missing_species,
          methode = methode,
          CS_version = CS_version,
          abstolint = abstolint,
          reltolint = reltolint
        )

    }
  }

  print_parameters_and_loglik(
    pars = pars,
    loglik = loglik,
    verbose = pars2[4],
    type = 'island_loglik'
  )

  return(loglik)

}

DAISIE_loglik_all_choosepar <- function(

  trparsopt,
  trparsfix,
  idparsopt,
  idparsfix,
  idparsnoshift,
  idparseq,
  pars2,
  datalist,
  methode,
  CS_version = 1,
  abstolint = 1E-16,
  reltolint = 1E-10

) {

  all_no_shift <- 6:10

  non_oceanic_option <- FALSE

  # So the 6th parameter is for non-oceanic!!
  if (

    max(idparsopt,-Inf) <= 6 &&
    max(idparsfix,-Inf) <= 6 &&
    (6 %in% idparsopt || 6 %in% idparsfix)

  ) {

    idparsnoshift <- 7:11
    all_no_shift <- 7:11
    non_oceanic_option <- TRUE

  }

  if (sum(idparsnoshift %in% (all_no_shift)) != 5) {

    # Type 2 stuff

    trpars1 <- rep(0, 11)

  } else {

    # No type 2 stuff

    trpars1 <- rep(0, 6)

    prop_type2_present <- which(idparsfix == 11)

    if(length(prop_type2_present) > 0) {

      trparsfix <- trparsfix[-prop_type2_present]
      idparsfix <- idparsfix[-prop_type2_present]

    }

  }

  trpars1[idparsopt] <- trparsopt

  if (length(idparsfix) != 0) {

    trpars1[idparsfix] <- trparsfix

  }

  # Again, if any type 2 stuff
  if (sum(idparsnoshift %in% all_no_shift) != 5) {

    trpars1[idparsnoshift] <- trpars1[idparsnoshift - 5]

  }

  if (max(trpars1) > 1 | min(trpars1) < 0) {

    # Those parameters should have been rescaled to max. 1
    loglik <- -Inf

  } else {

    # Rescale again
    pars1 <- trpars1 / (1 - trpars1)

    # If eequilibrium stuff
    if (pars2[6] > 0) {

      pars1 <- DAISIE_eq(datalist, pars1, pars2[-5])

      if (sum(idparsnoshift %in% all_no_shift) != 5) {

        pars1[idparsnoshift] <- pars1[idparsnoshift - 5]

      }
    }

    if (min(pars1) < 0 | (pars1[6] > 1 && non_oceanic_option == TRUE)) {

      # TODO

      loglik <- -Inf

    } else {

      # We are here
      loglik <- DAISIE::DAISIE_loglik_all(
        pars1 = pars1,
        pars2 = pars2,
        datalist = datalist,
        methode = methode,
        CS_version = CS_version,
        abstolint = abstolint,
        reltolint = reltolint
      )

    }

    if (is.nan(loglik) || is.na(loglik)) {

      message("There are parameter values used which cause numerical problems.")
      loglik <- -Inf

    }
  }

  return(loglik)

}

DAISIE_ML <- function(

  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  idparsnoshift = 6:10,
  res = 100,
  ddmodel = 0,
  cond = 0,
  eqmodel = 0,
  x_E = 0.95,
  x_I = 0.98,
  tol = c(1E-4, 1E-5, 1E-7),
  maxiter = 1000 * round((1.25) ^ length(idparsopt)),
  methode = "lsodes",
  optimmethod = "subplex",
  CS_version = 1,
  verbose = 0,
  tolint = c(1E-16, 1E-10),
  island_ontogeny = NA,
  jitter = 0,
  num_cycles = 1

) {

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

  # Set the parameters that must be at equilibrium depending on user choice
  if (eqmodel %in% c(1, 3, 13)) idparseq <- 2
  if (eqmodel %in% c(2, 4)) idparseq <-  4
  if (eqmodel %in% c(5, 15)) idparseq <- c(2, 4)

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

  # If...
  if (

    # ... no optimized parameter up to lambda_c2, and...
    max(idparsopt, -Inf) <= 6 &&

    # ... no fixed parameter up to lambda_c2, and...
    max(idparsfix, -Inf) <= 6 &&

    # lambda_c2 is either fixed or to be optimized
    (6 %in% idparsopt || 6 %in% idparsfix)

  ) {

    # Add probabilities of initial presence in the middle
    namepars <- c(namepars[1:5], "prob_init_pres", namepars[6:length(namepars)])

    # Add room for one more parameter
    max_idpars <- max_idpars + 1

    # Shift the indices of type 2 parameters
    all_no_shift <- all_no_shift + 1

    # Reset the indices of the parameters not-to-shift for type 2 species
    idparsnoshift <- all_no_shift

    # Prepare to add a new column of initial probabilities
    new_column <- NA
    names(new_column) <- "prob_init_pres"

    # Add a column of probabilities to the output, in the right position
    out2err <- add_column_to_dataframe(
      df = out2err,
      position = "lambda_a",
      column_to_insert = new_column
    )

  }

  # Tell the user what they are doing
  DAISIE:::print_ml_par_settings(

    namepars = namepars,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    idparsnoshift = idparsnoshift,
    all_no_shift = all_no_shift,
    verbose = TRUE

  )

  # idparsopt: those to optimize
  # idparsfix: those keep fixed
  # idparsnoshift: type 2 parameters not to consider different
  # idparseq: parameters at equilibrium

  # That should be all the possible indices no?
  idpars <- sort(c(idparsopt, idparsfix, idparsnoshift, idparseq))

  # Get the numbers of missing species for each clade
  missnumspec <- unlist(lapply(datalist, function(l) { l$missing_species }))

  # If the number of missing species for any clade is too large...
  if (max(missnumspec) > (res - 1)) {

    # Issue a warning
    warning(
      "The number of missing species is too large relative to the
       resolution of the ODE."
    )

    # And return an empty output
    return(out2err)

  }

  # Or if that number is a bit large (really?)
  if (max(missnumspec) > res / 10) {

    # Say it
    warning(
      "The number of missing species is quite low relative to the
        resolution of the ODE."
    )
  }

  # Check something
  if ((length(idpars) != max(idpars))) {

    warning("The parameters to be optimized and/or fixed are incoherent.")
    return(out2err)

  }


  # Some other check
  if (

    (!all(idpars == 1:max(idpars))) ||
    (length(initparsopt) != length(idparsopt)) ||
    (length(parsfix) != length(idparsfix))

  ) {

    warning("The parameters to be optimized and/or fixed are incoherent.")
    return(out2err)

  }

  if (length(idparseq) == 0) {

  } else {

    if (ddmodel == 3) {

      warning("Equilibrium optimization is not implemented for ddmodel = 3")

    } else {

      message(
        "You are assuming equilibrium. Extinction and/or immigration will
          be considered a function of the other parameters, the species
          pool size, the number of endemics,
          and/or the number of non-endemics"
      )

    }
  }

  # Rescale the initial guesses of optimised parameters (infinity becomes one)
  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1

  # Same for values of fixed parameters
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1

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

  # Compute the initial log-likelihood
  initloglik <- DAISIE:::DAISIE_loglik_all_choosepar(

    trparsopt = trparsopt,
    trparsfix = trparsfix,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    idparsnoshift = idparsnoshift,
    idparseq = idparseq,
    pars2 = pars2,
    datalist = datalist,
    methode = methode,
    CS_version = CS_version,
    abstolint = tolint[1],
    reltolint = tolint[2]

  )

  print_init_ll(initloglik = initloglik, verbose = verbose)

  if (initloglik == -Inf) {
    warning(
      "The initial parameter values have a likelihood that is equal to 0 or
       below machine precision. Try again with different initial values."
    )
    return(out2err)
  }

  out <- DDD::optimizer(
    optimmethod = optimmethod,
    optimpars = optimpars,
    fun = DAISIE_loglik_all_choosepar,
    trparsopt = trparsopt,
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
    jitter = jitter,
    num_cycles = num_cycles
  )

  if (out$conv != 0) {

    warning(
      "Optimization has not converged.
       Try again with different initial values."
    )

    out2 <- out2err
    out2$conv <- out$conv

    return(out2)

  }

  MLtrpars <- as.numeric(unlist(out$par))
  MLpars <- MLtrpars / (1 - MLtrpars)
  ML <- as.numeric(unlist(out$fvalues))

  if (sum(idparsnoshift %in% (all_no_shift)) != 5) {

    MLpars1 <- rep(0, 11)

  } else {

    MLpars1 <- rep(0, 6)

  }

  MLpars1[idparsopt] <- MLpars

  if (length(idparsfix) != 0) {

    MLpars1[idparsfix] <- parsfix

  }

  if (eqmodel > 0) {

    MLpars1 <- DAISIE_eq(datalist, MLpars1, pars2[-5])

  }

  if (MLpars1[3] > 10 ^ 7) {

    MLpars1[3] <- Inf

  }

  if (sum(idparsnoshift %in% (all_no_shift)) != 5) {

    if (length(idparsnoshift) != 0) {

      MLpars1[idparsnoshift] <- MLpars1[idparsnoshift - 5]

    }

    if (MLpars1[8] > 10 ^ 7) {

      MLpars1[8] <- Inf

    }

    out2 <- data.frame(

      lambda_c = MLpars1[1],
      mu = MLpars1[2],
      K = MLpars1[3],
      gamma = MLpars1[4],
      lambda_a = MLpars1[5],
      lambda_c2 = MLpars1[6],
      mu2 = MLpars1[7],
      K2 = MLpars1[8],
      gamma2 = MLpars1[9],
      lambda_a2 = MLpars1[10],
      prop_type2 = MLpars1[11],
      loglik = ML,
      df = length(initparsopt),
      conv = unlist(out$conv)

    )

    pars_to_print <- MLpars1[1:11]
    parnames <- c('lambda^c','mu','K','gamma','lambda^a','lambda^c2','mu2','K2','gamma2','lambda^a2','prop_type2')

  } else if (all(all_no_shift == 7:11)) {

    out2 <- data.frame(

      lambda_c = MLpars1[1],
      mu = MLpars1[2],
      K = MLpars1[3],
      gamma = MLpars1[4],
      lambda_a = MLpars1[5],
      prob_init_pres = MLpars1[6],
      loglik = ML,
      df = length(initparsopt),
      conv = unlist(out$conv)

    )

    pars_to_print <- MLpars1[1:6]
    parnames <- c("lambda^c", "mu", "K", "gamma", "lambda^a", "prob_init_pres")

  } else {

    out2 <- data.frame(

      lambda_c = MLpars1[1],
      mu = MLpars1[2],
      K = MLpars1[3],
      gamma = MLpars1[4],
      lambda_a = MLpars1[5],
      loglik = ML,
      df = length(initparsopt),
      conv = unlist(out$conv)

    )

    pars_to_print <- MLpars1[1:5]
    parnames <- c('lambda^c','mu','K','gamma','lambda^a')

  }

  print_parameters_and_loglik(
    pars = pars_to_print,
    loglik = ML,
    verbose = verbose,
    parnames = parnames,
    type = 'island_ML'
  )

  if (eqmodel > 0) {

    M <- calcMN(datalist, MLpars1)
    ExpEIN <- DAISIE_ExpEIN(datalist[[1]]$island_age, MLpars1, M) # nolint start

    message(
      paste0(
        "The expected number of endemics, non-endemics, and the total at ",
        "these parameters is: "
      ),
      paste(ExpEIN[[1]], ExpEIN[[2]], ExpEIN[[3]])
    ) # nolint end

  }

  return(invisible(out2))

}
