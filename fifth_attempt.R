rm(list = ls())

source("R/integrate_daisie.R")
source("R/right_hand_side.R")
source("R/get_nps.R")

pars <- list(lambda_c = 0.18, mu = 0.02, gamma = 0.02, lambda_a = 2)

integrate_daisie(island_age = -30, kpresent = 1L, pars, nmax = 10L, tcol = -10)

