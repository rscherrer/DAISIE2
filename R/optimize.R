# Function to find optimal parameter values for a function
optimize <- function(

  fun, pars, method = "simplex", tol = c(1e-04, 1e-04, 1e-06), maxiter = 1000,
  num_cycles = 1L, jitter = 0, verbose = FALSE, ...

) {

  # ...

  # TODO: Rename placeholder for function "fun" everywhere in the package.

  # Combine some parameters
  optimpars <- c(tol, maxiter)

  # Cap the number of cycles
  max_cycles <- ifelse(num_cycles == Inf, 10L, num_cycles)

  # TODO: Why would the number of cycles be infinite?

  cy <- 1

  fvalue <- rep(-Inf, max_cycles)

  out <- NULL

  while (cy <= max_cycles) {

    if (verbose & max_cycles > 1) message(paste("Cycle ", cy, "\n", sep = ""))

    if (method == "simplex") {

      # Note: control must contain reltolx, reltolf, abstolx, maxiter and verbose for simplex

      outnew <- suppressWarnings(with(control, simplex(

        fun, pars, reltolx, reltolf, abstolx, maxiter, ...

      )))

    } else if (method == "subplex") {

      # Note: control must contain reltol, maxit and jitter for subplex

      minfun1 <- function(fun, pars, ...) { return(-fun(pars = pars, ...)) }

      pars[pars == 0.5] <- 0.5 - control$jitter

      outnew <- suppressWarnings(subplex::subplex(

        par = pars, fn = minfun1, control = control, fun = fun, ...

      ))

      outnew <- list(par = outnew$par, fvalues = -outnew$value, conv = outnew$convergence)

    }

    if (cy > 1 & (any(is.na(outnew$par)) | any(is.nan(outnew$par)) |
                  is.na(outnew$fvalues) | is.nan(outnew$fvalues) |
                  outnew$conv != 0)) {

      if (verbose) message("The last cycle failed; second last cycle result is returned.\n")

      return(out)

    } else {

      out <- outnew
      pars <- out$par
      fvalue[cy] <- out$fvalues

    }

    if (cy > 1) {

      if (abs(fvalue[cy] - fvalue[cy - 1]) < optimpars[3]) {

        if (verbose & cy < max_cycles) message("No more cycles needed.\n")
        cy <- max_cycles

      } else if (cy == max_cycles) {

        message("More cycles in optimization recommended.\n")

      }

    }

    cy <- cy + 1

  }

  return(out)

}
