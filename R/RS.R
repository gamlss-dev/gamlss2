## Rigby and Stasinopoulos algorithm.
RS <- function(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
{
  ## Number of observations.
  n <- if(is.null(dim(y))) length(y) else nrow(y)

  ## Parameter names. FIXME: TRUE/FALSE?
  np <- family$names

  ## Initialize predictors.
  etastart <- initialize_eta(y, family, n, TRUE)

  ## Starting values.
  cstart <- NULL
  if(!missing(start)) {
    if(!inherits(start, "coef.gamlss2")) {
      if(!is.null(start)) {
        if(inherits(start, "list")) {
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
                etastart[[j]] <- rep(make.link2(family$links[j])$linkfun(start[[j]]), n)
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
  stop.eps <- eps
  eps <- eps + 1

  ## The step length control.
  if(is.null(control$step))
    control$step <- 1
  if((control$step > 1) | (control$step < 0))
    control$step <- 1

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
            fit[[j]]$fitted.values <- drop(x %*% fit[[j]]$coefficients)
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
    if(nes[[j]])
      etastart[[j]] <- eta[[j]]
  }

  ## Null deviance.
  dev0 <- -2 * family$loglik(y, family$map2par(eta))

  ## Estimate intercept only model first.
  if(isTRUE(control$nullmodel) & length(xterms)) {
    beta <- ieta <- list()
    for(j in np) {
      beta[[j]] <- as.numeric(fit[[j]]$coefficients["(Intercept)"])
      ieta[[j]] <- rep(beta[[j]], n)
    }
    beta <- unlist(beta)

    if(!any(is.na(beta))) {
      lli <- family$loglik(y, family$map2par(ieta))

      ridge <- control$ridge
      if(is.null(ridge))
        ridge <- 1e-05

      fn_ll <- function(par) {
        for(j in np)
          ieta[[j]] <- rep(par[j], n)
        ll <- family$loglik(y, family$map2par(ieta)) - ridge * sum(par^2)
        return(-ll)
      }

      opt <- try(nlminb(beta, objective = fn_ll), silent = TRUE)

      if(!inherits(opt, "try-error")) {
        if(-opt$objective > lli) {
          beta <- opt$par
          dev0 <- 2 * opt$objective
          if(isTRUE(control$initialize) & missing(start)) {
            for(j in np) {
              fit[[j]]$coefficients["(Intercept)"] <- beta[j]
              fit[[j]]$fitted.values <- drop(x[, "(Intercept)"] * fit[[j]]$coefficients["(Intercept)"])
              eta[[j]] <- fit[[j]]$fitted.values
            }
          }
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

  ## Track runtime.
  tstart <- proc.time()

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

  ## Start outer loop.
  while((eps[1L] > stop.eps[1L]) & iter[1L] < maxit[1L]) {
    ## Old log-likelihood.
    if(is.null(weights)) {
      llo0 <- family$loglik(y, family$map2par(eta))
    } else {
      llo0 <- sum(family$d(y, family$map2par(eta), log = TRUE) * weights, na.rm = TRUE)
    }

    ## Old predictors.
    if(iter[1L] >= CGk) {
      eta_old <- eta
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
      score <- deriv_checks(family$score[[j]](y, peta, id = j), is.weight = FALSE)
      hess <- deriv_checks(family$hess[[j]](y, peta, id = j), is.weight = TRUE)

      z <- if(iter[1L] > 0L) {
        eta[[j]] + 1 / hess * score
      } else {
        etastart[[j]] + 1 / hess * score
      }

      ## Start inner loop.
      while((eps[2L] > stop.eps[2L]) & iter[2L] < maxit[2L]) {
        ## Current log-likelihood.
        if(is.null(weights)) {
          ll0 <- family$loglik(y, family$map2par(eta))
        } else {
          ll0 <- sum(family$d(y, family$map2par(eta), log = TRUE) * weights, na.rm = TRUE)
        }
        ll02 <- ll0

        ## Cole and Green adjustment.
        if(iter[1L] >= CGk) {
          h <- grep(paste0(j, "."), names(family$hess), value = TRUE)
          if(length(h)) {
            ej <- sapply(strsplit(h, ".", fixed = TRUE), function(x) x[2L])
            adj <- 0.0
            for(l in seq_along(h)) {
              hess_l <- family$hess[[h[l]]](y, peta, id = j)  ## FIXME: deriv_checks()?
              adj <- adj + hess_l * (eta[[j]] - eta_old[[j]])
            }
            adj <- adj * hess
          }
        }

        ## Fit linear part.
        if(length(xterms[[j]])) {
          ## Compute partial residuals.
          eta[[j]] <- eta[[j]] - fit[[j]]$fitted.values
          if(iter[1L] >= CGk) {
            e <- (z - adj) - eta[[j]]
          } else {
            e <- z - eta[[j]]
          }

          ## Weights.
          wj <- if(is.null(weights)) hess else hess * weights

          ## Estimate weighted linear model.
          m <- lm.wfit(x[, xterms[[j]], drop = FALSE], e, wj, method = "qr")

          ## If linear model does not improve the fit, use ML.
          etai <- eta
          etai[[j]] <- etai[[j]] + m$fitted.values
          ll1 <- family$loglik(y, family$map2par(etai))

          if(ll1 < ll02) {
            ll <- function(par) {
              eta[[j]] <- eta[[j]] + drop(x[, xterms[[j]], drop = FALSE] %*% par)
              -family$loglik(y, family$map2par(eta)) + ridge * sum(par^2)
            }
            opt <- try(optim(coef(m), fn = ll, method = "BFGS"), silent = TRUE)
            if(!inherits(opt, "try-error")) {
              m$coefficients <- opt$par
              m$fitted.values <- drop(x[, xterms[[j]], drop = FALSE] %*% opt$par)
            }
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

          if(ll1 > ll02) {
            ## Update predictor.
            fit[[j]]$fitted.values <- m$fitted.values
            fit[[j]]$coefficients <- m$coefficients
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
            if(iter[1L] >= CGk) {
              e <- (z - adj) - eta[[j]]
            } else {
              e <- z - eta[[j]]
            }

            ## Additive model term fit.
            fs <- if(is.null(weights)) {
              special.wfit(specials[[k]], e, hess, y, eta, j, family, control,
                transfer = sfit[[j]][[k]]$transfer, iter = iter)
            } else {
              special.wfit(specials[[k]], e, hess * weights, y, eta, j, family, control,
                transfer = sfit[[j]][[k]]$transfer, iter = iter)
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

            if(ll1 > ll02) {
              ## Update predictor.
              sfit[[j]][[k]] <- fs
              sfit[[j]][[k]]$selected <- TRUE
              ll02 <- ll1
              ## sfit[[j]][[k]]$residuals <- z - etai[[j]] + fs$fitted.values ## FIXME: do we need this?
            } #else {
              #if(isTRUE(sfit[[j]][[k]]$selected)) {
                #sfit[[j]][[k]] <- fs
                #sfit[[j]][[k]]$selected <- TRUE
                #ll02 <- ll1
              #}
            #}

            eta[[j]] <- eta[[j]] + sfit[[j]][[k]]$fitted.values

            if(iter[1L] < 1L)
              etastart[[j]] <- eta[[j]]
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

    ## Warning if deviance is increasing?
    if((llo1 < llo0) & (iter[1L] > 0)) {
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
          if(!sfit[[j]][[i]]$selected)
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

  rval <- list(
    "fitted.values" = as.data.frame(eta),
    "fitted.specials" = sfit,
    "fitted.linear" = fit,
    "coefficients" = coef_lin,
    "elapsed" = elapsed, "iterations" = iter[1L],
    "logLik" = llo1, "control" = control,
    "nobs" = length(eta[[1L]]),
    "deviance" = -2 * llo1,
    "null.deviance" = dev0,
    "dev.reduction" = abs((dev0 - (-2 * llo1)) / dev0),
    "nullmodel" = control$nullmodel
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
  eta <- list()
  for(j in family$names)
    eta[[j]] <- rep(0.0, nobs)
  if(is.null(family$initialize) | !initialize)
    return(eta)
  for(j in family$names) {
    if(!is.null(family$initialize[[j]])) {
      linkfun <- make.link2(family$links[j])$linkfun
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

