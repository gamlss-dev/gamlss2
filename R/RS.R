RS <- function(x, y, specials, family, offsets, weights, xterms, sterms, control)
{
  ## Number of observations.
  n <- if(is.null(dim(y))) length(y) else nrow(y)

  ## Parameter names. FIXME: TRUE/FALSE?
  np <- family$names

  ## Initialize predictors.
  eta <- initialize_eta(y, family, n)

  ## Set control parameters.
  ## Stopping criterion.
  eps <- control$eps
  if(is.null(eps))
    eps <- sqrt(.Machine$double.eps)
  if(length(eps) < 2)
    eps <- c(eps, eps)
  stop.eps <- eps
  eps <- eps + 1

  ## Maximum number of backfitting iterations.
  maxit <- control$maxit
  if(is.null(maxit))
    maxit <- 400L
  if(length(maxit) < 2L)
    maxit <- c(maxit, 1L)

  ## Initialize fitted values for each model term.
  fit <- list()
  for(j in np) {
    fit[[j]] <- list()
    if(length(xterms[[j]]))
      fit[[j]]$linear <- list("fitted.values" = rep(0.0, n))
    if(length(sterms)) {
      if(length(sterms[[j]])) {
        for(i in sterms[[j]])
          fit[[j]][[i]] <- list("fitted.values" = rep(0.0, n))
      }
    }
  }

  ## Track runtime.
  tstart <- proc.time()

  ## Track iterations
  iter <- c(0, 0)

  ## For printing.
  if(control$flush) {
    control$flush <- interactive()
  }

  ## Start outer loop.
  while((eps[1L] > stop.eps[1L]) & iter[1L] < maxit[1L]) {
    ## Old log-likelihood.
    if(is.null(weights)) {
      llo0 <- family$loglik(y, family$map2par(eta))
    } else {
      llo0 <- sum(family$d(y, family$map2par(eta), log = TRUE) * weights, na.rm = TRUE)
    }

    for(j in np) {
      ## Outer loop working response and weights.
      peta <- family$map2par(eta)
      score <- family$score[[j]](y, peta, id = j)
      hess <- family$hess[[j]](y, peta, id = j)
      
      z <- eta[[j]] + 1 / hess * score

      ## Overwrite eta once.
      if(iter[1L] < 1) {
        eta[[j]] <- rep(0.0, n)
        if(nrow(offsets) > 0) {
          if(!is.null(offsets[[j]]))
             eta[[j]] <- eta[[j]] + offsets[[j]]
        }
      }

      ## Start inner loop.
      while((eps[2L] > stop.eps[2L]) & iter[2L] < maxit[2L]) {
        ## Current log-likelihood.
        if(is.null(weights)) {
          ll0 <- family$loglik(y, family$map2par(eta))
        } else {
          ll0 <- sum(family$d(y, family$map2par(eta), log = TRUE) * weights, na.rm = TRUE)
        }

        ## Fit linear part.
        if(length(xterms[[j]])) {
          ## Compute partial residuals.
          eta[[j]] <- eta[[j]] - fit[[j]]$linear$fitted.values
          e <- z - eta[[j]]

          ## Estimate weighted linear model.
          m <- if(is.null(weights)) {
            lm.wfit(x[, xterms[[j]], drop = FALSE], e, hess, method = "qr")
          } else {
            lm.wfit(x[, xterms[[j]], drop = FALSE], e, hess * weights, method = "qr")
          }

          ## Step length control.
          if(control$step != 1) {
            if((iter[1L] > 1L) | (iter[2L] > 1L))
              m$fitted.values <- control$step * m$fitted.values + (1 - control$step) * fit[[j]]$linear$fitted.values
          }

          ## Update predictor.
          eta[[j]] <- eta[[j]] + m$fitted.values
          fit[[j]]$linear$fitted.values <- m$fitted.values
          fit[[j]]$linear$coefficients <- m$coefficients
        }

        ## New log-likelihood.
        if(is.null(weights)) {
          ll1 <- family$loglik(y, family$map2par(eta))
        } else {
          ll1 <- sum(family$d(y, family$map2par(eta), log = TRUE) * weights, na.rm = TRUE)
        }

        ## Stopping criterion.
        eps[2L] <- abs((ll1 - ll0) / ll0)

        ## Update working response.
        if(eps[2L] > stop.eps[2L]) {
          peta <- family$map2par(eta)
          score <- family$score[[j]](y, peta, id = j)
          hess <- family$hess[[j]](y, peta, id = j)
          z <- eta[[j]] + 1 / hess * score
        }

        ## Update inner loop iterator.
        iter[2L] <- iter[2L] + 1L
      }

      ## Reset inner iterator and stopping criterion.
      iter[2L] <- 0
      eps[2L] <- stop.eps[2L] + 1
    }

    ## New log-likelihood.
    if(is.null(weights)) {
      llo1 <- family$loglik(y, family$map2par(eta))
    } else {
      llo1 <- sum(family$d(y, family$map2par(eta), log = TRUE) * weights, na.rm = TRUE)
    }

    ## Stopping criterion.
    eps[1L] <- abs((llo1 - llo0) / llo0)

    ## Update outer iterator.
    iter[1L] <- iter[1L] + 1L

    ## Print current state.
    if(control$trace) {
      if(iter[1L] > 1) {
        if(control$flush)
          cat("\r")
      }
      cat("GAMLSS-RS iteration ", iter[1L], ": Global Deviance = ",
        format(round(-2 * llo1, 4)), if(control$flush) NULL else "\n", sep = "")
    }
  }

  if(control$flush)
    cat("\n")

  ## Runtime.
  elapsed <- (proc.time() - tstart)["elapsed"]

  ## Extract linear/specials parts.
  coef_lin <- fit_specials <- list()
  for(j in np) {
    if(length(xterms[[j]])) {
      coef_lin[[j]] <- fit[[j]]$linear$coefficients
    }
  }

  rval <- list("fitted.values" = as.data.frame(eta),
    "coefficients" = coef_lin, "specials" = fit_specials,
    "elapsed" = elapsed, "iterations" = iter[1L],
    "logLik" = llo1)

  rval
}

## Function to initialize predictors.
initialize_eta <- function(y, family, nobs)
{
  eta <- list()
  for(j in family$names)
    eta[[j]] <- rep(0.0, nobs)
  if(is.null(family$initialize))
    return(eta)
  for(j in family$names) {
    if(!is.null(family$initialize[[j]])) {
      linkfun <- make.link2(family$links[j])$linkfun
      eta[[j]] <- linkfun(family$initialize[[j]](y))
      eta[[j]] <- rep(eta[[j]], length.out = nobs)
    }
  }
  return(eta)
}

