################################################################################
################################################################################
################################################################################
################################################################################
RS <- function(x, y, specials, family, offsets, weights, xterms, sterms, control)
{
  ## Number of observations.
  n <- if(is.null(dim(y))) length(y) else nrow(y)

  ## Parameter names. FIXME: TRUE/FALSE?
  np <- family$names

  ## Initialize predictors.
  eta <- initialize_eta(y, family, n)

  ## Check weights.
  if(!is.null(weights))
    weights <- as.numeric(weights)

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
  fit <- sfit <- list()
  for(j in np) {
    fit[[j]] <- list()
    if(length(xterms[[j]]))
      fit[[j]] <- list("fitted.values" = rep(0.0, n))
    if(length(sterms)) {
      if(length(sterms[[j]])) {
        sfit[[j]] <- list()
        for(i in sterms[[j]])
          sfit[[j]][[i]] <- list("fitted.values" = rep(0.0, n))
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
      score <- deriv_checks(family$score[[j]](y, peta, id = j), is.weight = FALSE)
      hess <- deriv_checks(family$hess[[j]](y, peta, id = j), is.weight = TRUE)

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
          eta[[j]] <- eta[[j]] - fit[[j]]$fitted.values
          e <- z - eta[[j]]

          ## Weights.
          wj <- if(is.null(weights)) hess else hess * weights

          ## Estimate weighted linear model.
          m <- lm.wfit(x[, xterms[[j]], drop = FALSE], e, wj, method = "qr")

          ## If linear model does not improve the fit, use ML.
          etai <- eta
          etai[[j]] <- etai[[j]] + m$fitted.values
          ll1 <- family$loglik(y, family$map2par(etai))
          if(ll1 < ll0) {
            ll <- function(par) {
              eta[[j]] <- eta[[j]] + drop(x[, xterms[[j]], drop = FALSE] %*% par)
              -family$loglik(y, family$map2par(eta))
            }
            start <- if(iter[1L] > 0 | iter[2L] > 0) coef(m) else rep(0, length(coef(m)))
            opt <- nlminb(start, objective = ll)
            m$coefficients <- opt$par
            m$fitted.values <- drop(x[, xterms[[j]], drop = FALSE] %*% opt$par)
          }
          
          ## Step length control.
          if(control$step < 1) {
            if(iter[1L] > 0 | iter[2L] > 0) {
              m$fitted.values <- control$step * m$fitted.values +
                (1 - control$step) * fit[[j]]$fitted.values
            }
          }

          etai <- eta
          etai[[j]] <- etai[[j]] + m$fitted.values
          ll1 <- family$loglik(y, family$map2par(etai))

          if(ll1 > ll0) {
            ## Update predictor.
            fit[[j]]$fitted.values <- m$fitted.values
            fit[[j]]$coefficients <- m$coefficients
            fit[[j]]$residuals <- z - etai[[j]] + m$fitted.values
          }
          eta[[j]] <- eta[[j]] + fit[[j]]$fitted.values
        }

        ## Fit specials part.
        if(length(sterms[[j]])) {
          for(k in sterms[[j]]) {
            ## Compute partial residuals.
            eta[[j]] <- eta[[j]] - sfit[[j]][[k]]$fitted.values
            e <- z - eta[[j]]

            ## Additive model term fit.
            fs <- if(is.null(weights)) {
              special.wfit(specials[[k]], e, hess, y, eta, j, family, control)
            } else {
              special.wfit(specials[[k]], e, hess * weights, y, eta, j, family, control)
            }

            ## Step length control.
            if(control$step < 1) {
              if(iter[1L] > 0 | iter[2L] > 0) {
                fs$fitted.values <- control$step * fs$fitted.values +
                  (1 - control$step) * sfit[[j]][[k]]$fitted.values
              }
            }

            etai <- eta
            etai[[j]] <- etai[[j]] + fs$fitted.values
            ll1 <- family$loglik(y, family$map2par(etai))

            if(ll1 > ll0) {
              ## Update predictor.
              sfit[[j]][[k]]$fitted.values <- fs$fitted.values
              sfit[[j]][[k]]$coefficients <- fs$coefficients
              sfit[[j]][[k]]$residuals <- z - etai[[j]] + fs$fitted.values
              sfit[[j]][[k]]$edf <- fs$edf
              sfit[[j]][[k]]$lambdas <- fs$lambdas
            }
            eta[[j]] <- eta[[j]] + sfit[[j]][[k]]$fitted.values
          }
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
          score <- deriv_checks(family$score[[j]](y, peta, id = j), is.weight = FALSE)
          hess <- deriv_checks(family$hess[[j]](y, peta, id = j), is.weight = TRUE)
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

    ## Warning if deviance is increasing.
    if((llo1 < llo0) & (iter[1L] > 0)) {
      warning("Deviance is increasing, maybe set argument step in gamlss2.control()!")
    }

    ## Update outer iterator.
    iter[1L] <- iter[1L] + 1L

    ## Print current state.
    if(control$trace) {
      if(iter[1L] > 1) {
        if(control$flush)
          cat("\r")
      }
      cat("GAMLSS-RS iteration ", fmt(iter[1L], nchar(as.character(maxit[1L])), digits = 0),
        ": Global Deviance = ", paste0(round(-2 * llo1, digits = 4), "   "),
        if(control$flush) NULL else "\n", sep = "")
    }
  }

  ## Runtime.
  elapsed <- (proc.time() - tstart)["elapsed"]

  if(control$trace & control$flush)
    cat("\n")

  ## Extract coefficients parts.
  coef_lin <- list()
  for(j in np) {
    if(length(xterms[[j]])) {
      coef_lin[[j]] <- fit[[j]]$coefficients
    }
  }

  rval <- list(
    "fitted.values" = as.data.frame(eta),
    "fitted.specials" = sfit,
    "fitted.linear" = fit,
    "coefficients" = coef_lin,
    "elapsed" = elapsed, "iterations" = iter[1L],
    "logLik" = llo1, "control" = control
  )

  class(rval) <- "gamlss2"

  rval
}
################################################################################
################################################################################
################################################################################
################################################################################
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

## Function to check values of score and hess vectors
deriv_checks <- function(x, is.weight = FALSE)
{
  x[is.na(x)] <- 1.490116e-08
  x[x > 1e+10] <- 1e+10
  if(is.weight) {
    x[x == 0] <- 1.490116e-08
    x[x < 0] <- -1 * x[x < 0]
    x[x < 1e-10] <- 1e-10
  } else {
    x[x < -1e+10] <- -1e+10
  }
  return(x)
}
################################################################################
################################################################################
################################################################################
################################################################################
## Formatting for printing.
fmt <- Vectorize(function(x, width = 8, digits = 2) {
  txt <- formatC(round(x, digits), format = "f", digits = digits, width = width)
  if(nchar(txt) > width) {
    txt <- strsplit(txt, "")[[1]]
    txt <- paste(txt[1:width], collapse = "", sep = "")
  }
  txt
})

fmt2 <- function(x, ...) {
  gsub(" ", "", fmt(x, ...))
}
################################################################################
################################################################################
################################################################################
################################################################################
