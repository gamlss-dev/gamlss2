## Rigby and Stasinopoulos algorithm.
RS <- function(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
{
  ## Number of observations.
  n <- if(is.null(dim(y))) length(y) else nrow(y)

  ## Parameter names. FIXME: TRUE/FALSE?
  np <- family$names

  ## Initialize predictors.
  etastart <- if(is.null(control$etastart)) TRUE else isTRUE(control$etastart)
  etastart <- initialize_eta(y, family, n, etastart)

  ## Starting values.
  cstart <- NULL
  if(missing(start))
    start <- NULL
  lp_start <- rep(FALSE, length(np))
  names(lp_start) <- np
  if(!is.null(start)) {
    if(!inherits(start, "coef.gamlss2")) {
      if(!is.null(start)) {
        if(inherits(start, c("gamlss2", "list"))) {
          if("fake_formula" %in% names(start)) {
            start <- fitted(start)
          } else {
            if(length(start[[1L]]) > 1L)
              start <- as.data.frame(start)
          }
        }
        if(inherits(start, c("data.frame", "matrix"))) {
          start <- as.data.frame(start)
          if(nrow(start) != n)
            stop("starting values have wrong number of observations!")
          for(j in np) {
            if(!is.null(start[[j]])) {
              lp_start[j] <- TRUE
              etastart[[j]] <- start[[j]]
            }
          }
        }
        if(inherits(start, c("list", "numeric"))) {
          start <- as.list(start)
          if(is.null(names(start)))
            names(start) <- rep(np, length.out = length(start))
          for(j in np) {
            if(!is.null(start[[j]])) {
              if(is.na(start[[j]]["(Intercept)"])) {
                etastart[[j]] <- rep(make.link2(family$links[[j]])$linkfun(start[[j]]), n)
              } else {
                etastart[[j]] <- rep(start[[j]], n)
              }
            }
          }
        }
      }
    } else {
      cstart <- unlist(start)
    }
  }

  ## Check weights.
  if(!is.null(weights))
    weights <- as.numeric(weights)

  ## Set control parameters.
  ## Stopping criterion.
  eps <- control$eps
  if(is.null(eps))
    eps <- 0.00001 ## sqrt(.Machine$double.eps)
  if(length(eps) < 2)
    eps <- c(eps, eps)
  eps <- rep(eps, length.out = 3)
  stop.eps <- eps
  eps <- eps + 1

  ## The step length control.
  if(is.null(control$step))
    control$step <- 1
  if((control$step > 1) | (control$step < 0))
    control$step <- 1
  if(is.null(control$autostep))
    control$autostep <- TRUE
  else
    control$autostep <- isTRUE(control$autostep)

  ## Maximum number of backfitting iterations.
  maxit <- control$maxit
  if(is.null(maxit))
    maxit <- 20L
  if(length(maxit) < 2L)
    maxit <- c(maxit, 20L)

  ## Fix some parameters?
  if(is.null(control$fixed)) {
    control$fixed <- rep(FALSE, length = length(np))
    names(control$fixed) <- np
  } else {
    if(is.null(names(control$fixed)))
      names(control$fixed) <- np[1:length(control$fixed)]
  }
  control$fixed <- as.list(control$fixed)
  for(j in np) {
    if(is.null(control$fixed[[j]]))
      control$fixed[[j]] <- FALSE
  }

  ## Ridge penalty?
  ridge <- isTRUE(control$ridge)
  penalty <- control$penalty
  if(is.null(penalty))
    penalty <- 1
  penalty <- rep(penalty, length.out = length(np))
  names(penalty) <- np

  ## Second fixed ridge penalty.
  lambda <- control$lambda
  if(is.null(lambda))
    lambda <- 1e-05

  ## Process offsets.
  if(!is.null(offsets)) {
    if(nrow(offsets) < 1L)
      offsets <- NULL
    else
      offsets <- as.data.frame(offsets)
  }
  
  ## Initialize fitted values for each model term.
  fit <- sfit <- eta <- nes <- list()
  for(j in np) {
    fit[[j]] <- list()
    eta[[j]] <- rep(0.0, n)
    nes[[j]] <- FALSE
    if(length(xterms[[j]])) {
      if("(Intercept)" %in% xterms[[j]]) {
        fit[[j]]$coefficients <- rep(0.0, length(xterms[[j]]))
        names(fit[[j]]$coefficients) <- xterms[[j]]
        fit[[j]]$coefficients["(Intercept)"] <- mean(etastart[[j]])
        fit[[j]]$fitted.values <- drop(x[, "(Intercept)"] * fit[[j]]$coefficients["(Intercept)"])
        if(!is.null(cstart)) {
          sj <- grep(paste0(j, ".p."), names(cstart), fixed = TRUE, value = TRUE)
          sj <- sj[sj %in% paste0(j, ".p.", xterms[[j]])]
          if(length(sj)) {
            fit[[j]]$coefficients[gsub(paste0(j, ".p."), "", sj)] <- as.numeric(cstart[sj])
            fit[[j]]$fitted.values <- drop(x[, names(fit[[j]]$coefficients), drop = FALSE] %*% fit[[j]]$coefficients)
            nes[[j]] <- TRUE
          }
        }
        eta[[j]] <- fit[[j]]$fitted.values
      } else {
        fit[[j]] <- list("fitted.values" = eta[[j]])
      }
    }
    if(length(sterms)) {
      if(length(sterms[[j]])) {
        sfit[[j]] <- list()
        for(i in sterms[[j]]) {
          sfit[[j]][[i]] <- list("fitted.values" = rep(0.0, n), "edf" = 0.0, "selected" = FALSE)
          if(!is.null(cstart)) {
            sj <- grep(paste0(j, ".s.", i), names(cstart), fixed = TRUE, value = TRUE)
            if(length(sj)) {
              if(!is.null(specials[[i]]$X)) {
                sfit[[j]][[i]]$fitted.values <- drop(specials[[i]]$X %*% cstart[sj])
                if(control$binning) {
                  sfit[[j]][[i]]$fitted.values <- sfit[[j]][[i]]$fitted.values[specials[[i]]$binning$match.index]
                }
                sfit[[j]][[i]]$selected <- TRUE
                eta[[j]] <- eta[[j]] + sfit[[j]][[i]]$fitted.values
                nes[[j]] <- TRUE
              }
            }
          }
        }
      }
    }
    if(!is.null(offsets)) {
      if(!is.null(offsets[[j]]))
        eta[[j]] <- eta[[j]] + offsets[[j]]
    }
    if(nes[[j]])
      etastart[[j]] <- eta[[j]]
  }

  ## Null deviance.
  dev0 <- -2 * family$logLik(y, family$map2par(etastart))

  ## Estimate intercept only model first.
  if(isTRUE(control$nullmodel) & length(unlist(xterms))) {
    nullmodel_ok <- TRUE
    beta <- ieta <- list()
    for(j in np) {
      beta[[j]] <- as.numeric(fit[[j]]$coefficients["(Intercept)"])
      if(length(beta[[j]]) < 1L)
        nullmodel_ok <- FALSE
      ieta[[j]] <- rep(beta[[j]], n)
      if(!is.null(offsets)) {
        if(!is.null(offsets[[j]]))
          ieta[[j]] <- ieta[[j]] + offsets[[j]]
      }
    }
    beta <- unlist(beta)

    if(!any(is.na(beta)) && nullmodel_ok) {
      lli <- family$logLik(y, family$map2par(ieta))

      fn_ll <- function(par) {
        for(j in np) {
          ieta[[j]] <- rep(par[j], n)
          if(!is.null(offsets)) {
            if(!is.null(offsets[[j]]))
              ieta[[j]] <- ieta[[j]] + offsets[[j]]
          }
        }
        ll <- family$logLik(y, family$map2par(ieta)) - lambda * sum(par^2)
        return(-ll)
      }

      opt <- try(nlminb(beta, objective = fn_ll), silent = TRUE)

      if(!inherits(opt, "try-error")) {
        if(-opt$objective > lli) {
          beta <- opt$par
          dev0 <- 2 * opt$objective
          if(isTRUE(control$initialize) & is.null(start)) {
            for(j in np) {
              fit[[j]]$coefficients["(Intercept)"] <- beta[j]
              fit[[j]]$fitted.values <- drop(x[, "(Intercept)"] * fit[[j]]$coefficients["(Intercept)"])
              eta[[j]] <- fit[[j]]$fitted.values
              if(!is.null(offsets)) {
                if(!is.null(offsets[[j]]))
                  eta[[j]] <- eta[[j]] + offsets[[j]]
              }
            }
          }
        } else {
          dev0 <- -2 * lli
        }
      }
    }
  }

  ## Use Cole and Green algorithm?
  CGk <- Inf
  if(!is.null(control$CG)) {
    if(!is.logical(control$CG))
      CGk <- as.integer(control$CG)
  }
  CG <- isTRUE(control$CG)
  if(CG)
    CGk <- 0L
  if(length(family$hess) < 2L)
    CGk <- Inf
  if(!any(grepl(".", names(family$hess), fixed = TRUE)))
    CGk <- Inf
  if(is.finite(CGk))
    eta_old <- eta
  if(length(maxit) < 3L) {
    if(is.finite(CGk))
      maxit <- c(maxit, 30)
    else
      maxit <- c(maxit, 3)
  }

  ## Track iterations
  iter <- c(0, 0)

  ## For printing.
  if(control$flush) {
    control$flush <- interactive()
  }

  if(control$trace) {
    if(!is.null(control$light)) {
      if(control$light)
        cat("Start estimation ...\n")
    }
  }

  step <- sapply(np, function(j) {
    rval <- list()
    if(length(xterms[[j]]))
      rval$xterms <- control$step
    if(length(sterms[[j]])) {
      rval$sterms <- rep(control$step, length(sterms[[j]]))
      names(rval$sterms) <- sterms[[j]]
    }
    rval
  }, simplify = FALSE)

  ## Start outer loop.
  while((eps[1L] > stop.eps[1L]) && (iter[1L] < maxit[1L])) {
    ## Old log-likelihood.
    if(is.null(weights)) {
      llo0 <- family$logLik(y, family$map2par(eta))
    } else {
      llo0 <- sum(family$pdf(y, family$map2par(eta), log = TRUE) * weights, na.rm = TRUE)
    }

    ## For CG.
    if(iter[1L] >= CGk) {
      eta_old <- if(iter[1L] > 0L) eta else etastart
      peta <- if(iter[1L] > 0L) {
        family$map2par(eta)
      } else {
        family$map2par(etastart)
      }
      zw_CG <- list()
      for(j in np) {
        zw_CG[[j]] <- z_weights(y, if(iter[1L] > 0L) eta[[j]] else etastart[[j]], peta, family, j)
      }
    }

    eps_outer <- 1
    iter_outer <- 0

    while((eps_outer > stop.eps[3L]) && (iter_outer < maxit[3L])) {
      if(iter[1L] >= CGk) {
        if(is.null(weights)) {
          outer_ll0 <- family$logLik(y, family$map2par(eta))
        } else {
          outer_ll0 <- sum(family$pdf(y, family$map2par(eta), log = TRUE) * weights, na.rm = TRUE)
        }
      }

      for(j in np) {
        ## Check if paramater is fixed.
        if(control$fixed[[j]])
          next

        ## Outer loop working response and weights.
        peta <- if(iter[1L] > 0L) {
          family$map2par(eta)
        } else {
          family$map2par(etastart)
        }

        ## Compute working response z and weights hess from family.
        ## Cole and Green adjustment.
        if(iter[1L] >= CGk) {
          h <- grep(paste0(j, "."), names(family$hess), value = TRUE)
          if(length(h)) {
            adj <- 0.0
            for(l in seq_along(h)) {
              parts <- strsplit(h[l], ".", fixed = TRUE)[[1]]
              k <- parts[2L]
              hess_l <- family$hess[[h[l]]](y, peta)
              if(!is.null(weights))
                hess_l <- hess_l * weights
              adj <- adj + hess_l * (eta[[k]] - eta_old[[k]])
            }
          }
          zw <- zw_CG[[j]]
          wj <- if(is.null(weights)) zw$weights else zw$weights * weights
          zw$z <- zw$z - adj / wj
        } else {
          zw <- z_weights(y, if(iter[1L] > 0L) eta[[j]] else etastart[[j]], peta, family, j)
        }

        ## Start inner loop.
        while((eps[2L] > stop.eps[2L]) && (iter[2L] < maxit[2L])) {
          ## Current log-likelihood.
          if(is.null(weights)) {
            ll0 <- family$logLik(y, family$map2par(eta))
          } else {
            ll0 <- sum(family$pdf(y, family$map2par(eta), log = TRUE) * weights, na.rm = TRUE)
          }
          ll02 <- ll0

          ## Fit linear part.
          if(length(xterms[[j]])) {
            ## Compute partial residuals.
            eta[[j]] <- eta[[j]] - fit[[j]]$fitted.values
            e <- zw$z - eta[[j]]

            ## Weights.
            wj <- if(is.null(weights)) zw$weights else zw$weights * weights

            ## Design matrix.
            Xj <- x[, xterms[[j]], drop = FALSE]

            ## Estimate weighted linear model.
            if(ridge) {
              m <- ridge.lm.wfit(Xj, e, wj, penalty = penalty[j], control)
              penalty[j] <- m$penalty
            } else {
              m <- lm.wfit(Xj, e, wj, method = "qr")
            }

            ## If linear model does not improve the fit, use ML.
            etai <- eta
            etai[[j]] <- etai[[j]] + m$fitted.values

            if(is.null(weights)) {
              ll1 <- family$logLik(y, family$map2par(etai))
            } else {
              ll1 <- sum(family$pdf(y, family$map2par(etai), log = TRUE) * weights, na.rm = TRUE)
            }

            if(ll1 < ll02 && isTRUE(control$backup)) {
              ll <- function(par) {
                etai <- eta
                etai[[j]] <- etai[[j]] + drop(Xj %*% par)
                -family$logLik(y, family$map2par(etai)) + lambda * sum(par^2)
              }
              warn <- getOption("warn")
              options("warn" = -1)
              opt <- try(optim(coef(m), fn = ll, method = "BFGS"), silent = TRUE)
              opt2 <- try(nlminb(coef(m), ll), silent = TRUE)
              options("warn" = warn)
              if(!inherits(opt2, "try-error")) {
                if(!inherits(opt, "try-error")) {
                  if(opt2$objective < opt$value) {
                    opt$par <- opt2$par
                  }
                } else {
                  opt <- opt2
                }
              }
              if(!inherits(opt, "try-error")) {
                m$coefficients <- opt$par
                m$fitted.values <- drop(Xj %*% opt$par)
                etai[[j]] <- etai[[j]] + m$fitted.values
                if(is.null(weights)) {
                  ll1 <- family$logLik(y, family$map2par(etai))
                } else {
                  ll1 <- sum(family$pdf(y, family$map2par(etai), log = TRUE) * weights, na.rm = TRUE)
                }
              }
            }
          
            ## Step length control.
            if((step[[j]]$xterms < 1 || control$autostep) && (ll1 < ll02)) {
              if(iter[1L] > 0 | iter[2L] > 0) {
                if(control$autostep) {
                  stepfun <- function(nu) {
                    b <- nu * m$coefficients + (1 - nu) * fit[[j]]$coefficients
                    f <- drop(Xj %*% b)
                    etai <- eta
                    etai[[j]] <- etai[[j]] + f
                    -family$logLik(y, family$map2par(etai))
                  }
                  s <- try(optimize(stepfun, lower = -1, upper = 1, tol = .Machine$double.eps^0.5), silent = TRUE)
                  if(-s$objective > ll02) {
                    step[[j]]$xterms <- s$minimum
                  }
                }
                m$coefficients <- step[[j]]$xterms * m$coefficients +
                  (1- step[[j]]$xterms) * fit[[j]]$coefficients
                m$fitted.values <- drop(Xj %*% m$coefficients)
              }
            }

            etai <- eta
            etai[[j]] <- etai[[j]] + m$fitted.values

            if(is.null(weights)) {
              ll1 <- family$logLik(y, family$map2par(etai))
            } else {
              ll1 <- sum(family$pdf(y, family$map2par(etai), log = TRUE) * weights, na.rm = TRUE)
            }

            if(ll1 > ll02) {
              ## Update predictor.
              fit[[j]]$fitted.values <- m$fitted.values
              fit[[j]]$coefficients <- m$coefficients
              if(!is.null(m$edf))
                fit[[j]]$edf <- m$edf
              XWX <- crossprod(Xj * sqrt(wj)) + diag(1e-6, ncol(Xj))
              fit[[j]]$vcov <- chol2inv(chol(XWX))
              colnames(fit[[j]]$vcov) <- rownames(fit[[j]]$vcov) <- colnames(Xj)
              ll02 <- ll1
            ## fit[[j]]$residuals <- z - etai[[j]] + m$fitted.values ## FIXME: do we need this?
            }

            eta[[j]] <- eta[[j]] + fit[[j]]$fitted.values

            if(iter[1L] < 1L)
              etastart[[j]] <- eta[[j]]
          }

          ## Fit specials part.
          if(length(sterms[[j]])) {
            for(k in sterms[[j]]) {
              ## Compute partial residuals.
              eta[[j]] <- eta[[j]] - sfit[[j]][[k]]$fitted.values
              e <- zw$z - eta[[j]]

              ## Additive model term fit.
              fs <- if(is.null(weights)) {
                special.wfit(specials[[k]], e, zw$weights, y, eta, j, family, control,
                  transfer = sfit[[j]][[k]]$transfer, iter = iter)
              } else {
                special.wfit(specials[[k]], e, zw$weights * weights, y, eta, j, family, control,
                  transfer = sfit[[j]][[k]]$transfer, iter = iter)
              }

              ## Step length control.
              if(step[[j]]$sterms[k] < 1) {
                if(iter[1L] > 0 | iter[2L] > 0) {
                  if(inherits(specials[[k]], c("mgcv.smooth", "X %*% b"))) {
                    fs$coefficients <- step[[j]]$sterms[k] * fs$coefficients +
                      (1 - step[[j]]$sterms[k]) * if(is.null(sfit[[j]][[k]]$coefficients)) 0 else sfit[[j]][[k]]$coefficients
                    fs$fitted.values <- drop(specials[[k]]$X %*% fs$coefficients)
                    if(control$binning)
                      fs$fitted.values <- fs$fitted.values[specials[[k]]$binning$match.index]
                  } else {
                    fs$fitted.values <- step[[j]]$sterms[k] * fs$fitted.values +
                      (1 - step[[j]]$sterms[k]) * sfit[[j]][[k]]$fitted.values
                  }
                }
              }

              etai <- eta
              etai[[j]] <- etai[[j]] + fs$fitted.values

              if(is.null(weights)) {
                ll1 <- family$logLik(y, family$map2par(etai))
              } else {
                ll1 <- sum(family$pdf(y, family$map2par(etai), log = TRUE) * weights, na.rm = TRUE)
              }

              if(ll1 > ll02) {
                ## Update predictor.
                sfit[[j]][[k]] <- fs
                sfit[[j]][[k]]$selected <- TRUE
                ll02 <- ll1
                ## sfit[[j]][[k]]$residuals <- z - etai[[j]] + fs$fitted.values ## FIXME: do we need this?
              } else {
                if(control$autostep) {
                  step[[j]]$sterms[k] <- step[[j]]$sterms[k] * 0.5
                }
              }

              eta[[j]] <- eta[[j]] + sfit[[j]][[k]]$fitted.values

              if(iter[1L] < 1L)
                etastart[[j]] <- eta[[j]]
            }
          }

          ## New log-likelihood.
          if(is.null(weights)) {
            ll1 <- family$logLik(y, family$map2par(eta))
          } else {
            ll1 <- sum(family$pdf(y, family$map2par(eta), log = TRUE) * weights, na.rm = TRUE)
          }

          ## Stopping criterion.
          eps[2L] <- abs(ll1 - ll0) / (abs(ll0) + 1e-08)

          ## Update working response.
          if((eps[2L] > stop.eps[2L]) && (iter[1L] < CGk)) {
            peta <- family$map2par(eta)
            zw <- z_weights(y, if(iter[1L] > 0L) eta[[j]] else etastart[[j]], peta, family, j)
          }

          ## Update inner loop iterator.
          iter[2L] <- iter[2L] + 1L
        }

        ## Reset inner iterator and stopping criterion.
        iter[2L] <- 0
        eps[2L] <- stop.eps[2L] + 1
      }

      ## For Cole and Green.
      iter_outer <- iter_outer + 1L
      if(iter[1L] >= CGk)
        eps_outer <- abs((ll1 - outer_ll0) / ll1)
    }

    ## New log-likelihood.
    if(is.null(weights)) {
      llo1 <- family$logLik(y, family$map2par(eta))
    } else {
      llo1 <- sum(family$pdf(y, family$map2par(eta), log = TRUE) * weights, na.rm = TRUE)
    }

    ## Stopping criterion.
    eps[1L] <- abs((llo1 - llo0) / llo0)

    ## Warning if deviance is increasing?
    if((llo1 < llo0) && (iter[1L] > 0)) {
      warning("Deviance is increasing, maybe set argument step in gamlss2_control()!")
    }

    ## Update outer iterator.
    iter[1L] <- iter[1L] + 1L

    ## Print current state.
    if(control$trace) {
      if(iter[1L] > 1) {
        if(control$flush) {
          cat('\r')
        }
      }
      itxt <- paste0(paste0("GAMLSS-", if(iter[1L] >= CGk) "CG" else "RS", " iteration "),
        fmt(iter[1L], nchar(as.character(maxit[1L])), digits = 0),
        ": Global Deviance = ", round(-2 * llo1, digits = 4),
        " eps = ", fmt(eps[1L], width = 8, digits = 8), "    ")
      cat(itxt, if(control$flush) NULL else "\n")
    }
  }

  if(control$trace & control$flush)
    cat("\n")

  ## Extract coefficients parts.
  coef_lin <- list()
  for(j in np) {
    if(length(xterms[[j]])) {
      coef_lin[[j]] <- fit[[j]]$coefficients
      ##coef_lin[[j]][is.na(coef_lin[[j]])] <- 0.0
    }
  }

  ## Check if special terms are never updated
  ## and remove fitted values if light = TRUE.
  if(length(sfit)) {
    dropj <- NULL
    for(j in names(sfit)) {
      if(length(sfit[[j]])) {
        drop <- NULL
        for(i in names(sfit[[j]])) {
          if(control$light) {
            sfit[[j]][[i]]$fitted.values <- NULL
          }
          if(!isTRUE(sfit[[j]][[i]]$selected))
            drop <- c(drop, i)
        }
        if(length(drop)) {
          sfit[[j]][drop] <- NULL
        }
        if(length(sfit[[j]]) < 1L)
          dropj <- c(dropj, j)
      }
    }
    if(length(dropj))
      sfit[dropj] <- NULL
  }

  if(ridge) {
    attr(fit, "edf") <- unlist(sapply(fit, function(x) x$edf))
  }

  ## Message if not converged due to NAs or Inf!
  d <- family$pdf(y, family$map2par(eta), log = TRUE)
  if(!is.null(weights))
    d <- d * weights
  if(any(is.na(d))) {
    warning("NA log-density values in the last iteration of the RS algorithm!")
  }
  if(any(!is.finite(d))) {
    warning("non-finite log-density values in the last iteration of the RS algorithm!")
  }

  rval <- list(
    "fitted.values" = as.data.frame(eta),
    "fitted.specials" = sfit,
    "fitted.linear" = fit,
    "coefficients" = coef_lin,
    "iterations" = iter[1L],
    "logLik" = llo1, "control" = control,
    "nobs" = length(eta[[1L]]),
    "deviance" = -2 * llo1,
    "null.deviance" = dev0,
    "dev.reduction" = (dev0 - (-2 * llo1)) / dev0,
    "nullmodel" = control$nullmodel,
    "stepsize" = step
  )

  class(rval) <- "gamlss2"

  rval
}

## Cole and Green flavor.
CG <- function(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
{
  control$CG <- TRUE
  RS(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
}

## Function to initialize predictors.
initialize_eta <- function(y, family, nobs, initialize)
{
  if(is.null(initialize)) {
    initialize <- TRUE
  } else {
    initialize <- isTRUE(initialize)
  }
  eta <- list()
  for(j in family$names)
    eta[[j]] <- rep(0.0, nobs)
  if(is.null(family$initialize) | !initialize)
    return(eta)
  for(j in family$names) {
    if(!is.null(family$initialize[[j]])) {
      linkfun <- make.link2(family$links[[j]])$linkfun
      eta[[j]] <- try(linkfun(family$initialize[[j]](y)), silent = TRUE)
      if(inherits(eta[[j]], "try-error")) {
        if(is.null(dim(y)))
          eta[[j]] <- linkfun(family$initialize[[j]](matrix(y, ncol = 1)))
      }
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
    x[(x == 0) | !is.finite(x)] <- 1.490116e-08
    x[x < 0] <- -1 * x[x < 0]
    x[x < 1e-10] <- 1e-10
  } else {
    x[x < -1e+10] <- -1e+10
  }
  return(x)
}

## Compute working response z and weights hess from family.
z_weights <- function(y, eta, peta, family, j)
{
  if(is.null(family$z_weights)) {
    score <- deriv_checks(family$score[[j]](y, peta, id = j), is.weight = FALSE)
    hess <- deriv_checks(family$hess[[j]](y, peta, id = j), is.weight = TRUE)
    z <- eta + 1 / hess * score
    return(list("z" = z, "weights" = hess))
  } else {
    return(family$z_weights(y, eta, peta, j))
  }
}

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

## Ridge regression linear model.
ridge.lm.wfit <- function(x, y, w, penalty, control)
{
  if(is.null(control$criterion))
    control$criterion <- "gaic"

  K <- control$K
  if(is.null(K))
    K <- 2.0

  ## Precompute weighted crossproducts using sqrt(w) once.
  sw <- sqrt(w)
  XW <- x * sw
  XWX <- crossprod(XW)
  XWy <- crossprod(XW, y * sw)

  nc <- ncol(x)
  i <- which(colnames(x) == "(Intercept)")
  I <- rep(1, nc)
  if(length(i))
    I[i] <- 0.0

  only_itcpt <- all(colnames(x) == "(Intercept)")

  n <- length(w)

  ## Function to find optimum ridge penalty.
  fp <- function(pen, rf = FALSE) {

    ## Penalty matrix.
    S <- diag(I * pen)
    if(!length(S))
      S <- 0.0

    ## Precision matrix Q = X'WX + S.
    Q <- XWX + S
    Q <- Q + diag(1e-08, ncol(Q))

    ## Cholesky factorization.
    cholQ <- try(chol(Q), silent = TRUE)
    if(inherits(cholQ, "try-error")) {
      Q <- Q + diag(1e-05, ncol(Q))
      cholQ <- chol(Q)
    }

    ## b = Q^{-1} X'Wy.
    b <- backsolve(cholQ, forwardsolve(t(cholQ), XWy))
    b <- drop(b)

    fit <- drop(x %*% b)

    ## EDF = tr(X'WX Q^{-1}).
    Tmat <- backsolve(cholQ, forwardsolve(t(cholQ), XWX))
    edf <- sum(diag(Tmat))

    ## Guard: edf can get numerically >= n in extreme cases.
    if(!is.finite(edf))
      edf <- nc
    if(edf > (n - 1e-08))
      edf <- n - 1e-08

    if(rf) {
      names(b) <- colnames(x)
      return(list("coefficients" = b, "fitted.values" = fit, "edf" = edf,
        "penalty" = pen, "vcov" = chol2inv(cholQ), "df" = n - edf))
    } else {
      rss <- sum(w * (y - fit)^2)

      rval <- switch(tolower(control$criterion),
        "gcv"  = rss * n / (n - edf)^2,
        "aic"  = rss + 2 * edf,
        "gaic" = rss + K * edf,
        "aicc" = rss + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1),
        "bic"  = rss + log(n) * edf
      )

      return(rval)
    }
  }

  if(!only_itcpt) {
    ## Optimize ridge penalty over a modest bracket.
    opt <- nlminb(penalty, objective = fp, lower = penalty / 10, upper = penalty * 10)
  } else {
    opt <- list(par = 0.0)
  }

  rval <- fp(opt$par, rf = TRUE)

  return(rval)
}

