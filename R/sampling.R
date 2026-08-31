## Function to compute multivariate normal samples
## given the maximum likelihood estimator.
sampling <- function(object, R = 100, antithetic = TRUE, ...)
{
  V <- vcov(object, ...)

  Cv <- tryCatch(
    chol(V),
    error = function(e) NULL
  )

  if(is.null(Cv)) {
    V <- as.matrix(Matrix::nearPD(V)$mat)
    Cv <- chol(V)
  }

  cb <- coef(object, dropall = FALSE, ...)
  p <- length(cb)

  ## Antithetic sampling gives exact symmetry around zero
  ## if R is even. If R is odd, add one extra draw and drop later.
  if(antithetic) {
    R0 <- ceiling(R / 2)
    Z <- matrix(rnorm(R0 * p), nrow = p, ncol = R0)
    Z <- cbind(Z, -Z)
    Z <- Z[, seq_len(R), drop = FALSE]
  } else {
    Z <- matrix(rnorm(R * p), nrow = p, ncol = R)
  }

  ## Transform standard normal draws.
  sc <- t(Cv) %*% Z

  ## Force exact zero column means coefficient-wise.
  sc <- sweep(sc, 1L, rowMeans(sc), "-")

  ## Add coefficient vector.
  sc <- sweep(sc, 1L, cb, "+")

  sc <- t(sc)
  colnames(sc) <- names(cb)

  sc
}

## Bayesian GAMLSS sampler function. Unbounded slice sampler!?
BS <- function(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
{
  ## Number of observations.
  n <- if(is.null(dim(y))) length(y) else nrow(y)

  ## Parameter names. FIXME: TRUE/FALSE?
  np <- family$names

  ## Initialize predictors.
  etastart <- initialize_eta(y, family, n, TRUE)

  ## Number of iterations.
  n.iter <- control$n.iter
  if(is.null(n.iter))
    n.iter <- 1200L

  ## Burnin period.
  burnin <- control$burnin
  if(is.null(burnin))
    burnin <- 200L

  ## Thinning samples.
  thin <- control$thin
  if(is.null(thin))
    thin <- 1L

  ## Type conversion.
  n.iter   <- as.integer(n.iter)
  burnin   <- as.integer(burnin)
  thin <- as.integer(thin)

  ## Basic sanity checks.
  if(is.na(n.iter) || n.iter <= 0L)
    stop("n.iter must be a positive integer.")

  if(is.na(burnin) || burnin < 0L)
    stop("burnin must be a non-negative integer.")

  if(is.na(thin) || thin <= 0L)
    stop("thin must be a positive integer.")

  ## Logical consistency (adaptive).
  if(burnin >= n.iter)
    burnin <- 0L

  ## Number of saved iterations (adaptive).
  nsave <- (n.iter - burnin) %/% thin

  if(nsave <= 0L) {
    thin <- 1L
    nsave <- n.iter - burnin
  }

  if(nsave <= 0L) {
    burnin <- 0L
    thin <- 1L
    nsave <- n.iter
  }

  ## Numbers of samples to save.
  iterthin <- seq.int(burnin + 1L, n.iter, by = thin)
  nsave <- length(iterthin)

  ## Starting values [same in RS()].
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
  fit <- sfit <- eta <- nes <- samples <- list()
  for(j in np) {
    samples[[j]] <- list()
    fit[[j]] <- list()
    eta[[j]] <- rep(0.0, n)
    nes[[j]] <- FALSE
    if(length(xterms[[j]])) {
      samples[[j]]$p <- matrix(NA, nrow = nsave, ncol = length(xterms[[j]]) + 1L)
      colnames(samples[[j]]$p) <- c(xterms[[j]], "alpha")
      fit[[j]]$coefficients <- rep(0.0, length(xterms[[j]]))
      names(fit[[j]]$coefficients) <- xterms[[j]]
      if("(Intercept)" %in% xterms[[j]]) {
        fit[[j]]$coefficients["(Intercept)"] <- mean(etastart[[j]])
        fit[[j]]$fitted.values <- drop(x[, "(Intercept)"] * fit[[j]]$coefficients["(Intercept)"])
      } else {
        fit[[j]]$fitted.values <- eta[[j]]
      }
      if(!is.null(cstart)) {
        sj <- grep(paste0(j, ".p."), names(cstart), fixed = TRUE, value = TRUE)
        sj <- sj[sj %in% paste0(j, ".p.", xterms[[j]])]
        if(length(sj)) {
          fit[[j]]$coefficients[gsub(paste0(j, ".p."), "", sj)] <- as.numeric(cstart[sj])
          fit[[j]]$fitted.values <- drop(x[, xterms[[j]], drop = FALSE] %*% fit[[j]]$coefficients)
          nes[[j]] <- TRUE
        }
      }
      eta[[j]] <- fit[[j]]$fitted.values
    }
    if(length(sterms)) {
      if(length(sterms[[j]])) {
        samples[[j]]$s <- list()
        sfit[[j]] <- list()
        for(i in sterms[[j]]) {
          if(!inherits(specials[[i]], c("mgcv.smooth", "mcmc"))) {
            stop("only mgcv and mcmc smooth terms are allowed!")
          }
          sfit[[j]][[i]] <- list(
            "fitted.values" = rep(0.0, n),
            "edf" = 0.0,
            "coefficients" = setNames(rep(0.0, ncol(specials[[i]]$X)),
              paste0(j, ".s.", i, ".", seq_len(ncol(specials[[i]]$X)))),
            "tau" = setNames(rep(0.001, length(specials[[i]]$S)),
              paste0(j, ".s.", i, ".tau", seq_along(specials[[i]]$S)))
          )
          samples[[j]]$s[[i]] <- matrix(NA, nrow = nsave,
            ncol = ncol(specials[[i]]$X) + length(specials[[i]]$S) + 2L)
          if(!is.null(cstart)) {
            prefix <- paste0(j, ".s.", i, ".")
            sj <- names(cstart)[startsWith(names(cstart), prefix)]
            sjb <- sj[!grepl(".lambda", sj, fixed = TRUE)]
            sjb <- sjb[!grepl(".tau", sjb, fixed = TRUE)]
            if(length(sjb)) {
              sfit[[j]][[i]]$fitted.values <- drop(specials[[i]]$X %*% cstart[sjb])
              sfit[[j]][[i]]$coefficients <- cstart[sjb]
              if(control$binning) {
                sfit[[j]][[i]]$fitted.values <- sfit[[j]][[i]]$fitted.values[specials[[i]]$binning$match.index]
              }
              sfit[[j]][[i]]$selected <- TRUE
              eta[[j]] <- eta[[j]] + sfit[[j]][[i]]$fitted.values
              nes[[j]] <- TRUE
            }
            sjl <- sj[grepl(".lambda", sj, fixed = TRUE)]
            if(length(sjl)) {
              sfit[[j]][[i]]$tau <- 1 / cstart[sjl]
              names(sfit[[j]][[i]]$tau) <- gsub("lambda", "tau", names(sfit[[j]][[i]]$tau))
            } else {
              sjt <- sj[grepl(".tau", sj, fixed = TRUE)]
              if(length(sjt)) {
                sfit[[j]][[i]]$tau <- cstart[sjt]
                names(sfit[[j]][[i]]$tau) <- gsub("lambda", "tau", names(sfit[[j]][[i]]$tau))
              }
            }
          }

          ## Assign a prior.
          if(is.null(specials[[i]]$prior))
            specials[[i]]$prior <- prior(specials[[i]])
        }
      }
    }
    if(nes[[j]])
      etastart[[j]] <- eta[[j]]
  }

  if(!is.null(control$fixed)) {
    for(j in np) {
      if(control$fixed[[j]]) {
        link <- make.link2(family$links[[j]])
        fit[[j]]$coefficients["(Intercept)"] <- link$linkfun(control$fixed[[j]])
        eta[[j]] <- rep(fit[[j]]$coefficients["(Intercept)"], n)
        fit[[j]]$fitted.values <- eta[[j]]
        etastart[[j]] <- eta[[j]]
      }
    }
  }

  ## Null deviance.
  dev0 <- -2 * family$log_likelihood(par = family$map2par(etastart), y = y)

  ## Estimate intercept only model first.
  if(isTRUE(control$nullmodel) & length(xterms)) {
    beta <- ieta <- list()
    for(j in np) {
      beta[[j]] <- as.numeric(fit[[j]]$coefficients["(Intercept)"])
      ieta[[j]] <- rep(beta[[j]], n)
    }
    beta <- unlist(beta)

    if(!any(is.na(beta))) {
      lli <- family$log_likelihood(par = family$map2par(ieta), y = y)

      fn_ll <- function(par) {
        for(j in np) {
          if(control$fixed[[j]])
            par[j] <- beta[j]
          ieta[[j]] <- rep(par[j], n)
        }
        ll <- family$log_likelihood(par = family$map2par(ieta), y = y) - 1e-05 * sum(par^2)
        return(-ll)
      }

      opt <- try(nlminb(beta, objective = fn_ll), silent = TRUE)

      if(!inherits(opt, "try-error")) {
        if(-opt$objective > lli) {
          beta <- opt$par
          dev0 <- 2 * opt$objective
          if(isTRUE(control$initialize) & (missing(start) | is.null(start))) {
            for(j in np) {
              fit[[j]]$coefficients["(Intercept)"] <- beta[j]
              fit[[j]]$fitted.values <- drop(x[, "(Intercept)"] * fit[[j]]$coefficients["(Intercept)"])
              eta[[j]] <- fit[[j]]$fitted.values
            }
          }
        }
        ## Null deviance.
        dev0 <- -2 * family$log_likelihood(par = family$map2par(eta), y = y)
      }
    }
  }

  ## Start MCMC.
  if(control$trace) {
    if(!is.null(control$light)) {
      if(control$light)
        cat("Start sampling ...\n")
    }
  }

  ## Priors.
  priors <- control$priors

  if(is.null(priors)) {
    priors$p <- function(parameters) {
      sum(dnorm(parameters, sd = 1000, log = TRUE))
    }
  }

  ## Tracking.
  track <- list()
  track$logLik <- rep(NA_real_, nsave)
  track$deviance <- rep(NA_real_, nsave)
  track$eta <- list()
  for(j in np) {
    track$eta[[j]] <- rep(0.0, n)
  }

  ## Start time etc.
  ptm <- proc.time()
  step <- 20
  nstep <- step
  step <- floor(n.iter / step)
  isave <- 1L

  for(iter in seq_len(n.iter)) {
    do_save <- (isave <= length(iterthin) && iter == iterthin[isave])
    for(j in np) {
      ## Check if paramater is fixed.
      if(control$fixed[[j]]) {
        if(do_save) {
          samples[[j]]$p[isave, ] <- c(fit[[j]]$coefficients, 1)
        }
        next
      }

      ## Sampling linear part.
      if(length(xterms[[j]])) {
        ## Get parameters.
        peta <- family$map2par(eta)

        ## Compute old log-likelihood.
        pibeta <- family$log_likelihood(par = peta, y = y)

        ## Old parameters.
        b0 <- fit[[j]]$coefficients

        ## Log-prior.
        p1 <- priors$p(b0)

        ## Derivatives.
        score <- deriv_checks(family$score[[j]](par = peta, y = y, id = j), is.weight = FALSE)
        hessian <- deriv_checks(family$hessian[[j]](par = peta, y = y, id = j), is.weight = TRUE)

        ## Working response.
        z <- eta[[j]] + 1 / hessian * score

        ## Compute residuals.
        eta2 <- eta[[j]] <- eta[[j]] - fit[[j]]$fitted.values
        e <- z - eta2

        ## Weights.
        wj <- if(is.null(weights)) hessian else hessian * weights

        ## Compute mean and precision.
        Xj <- x[, xterms[[j]], drop = FALSE]
        XW <- Xj * sqrt(wj)
        XWX <- crossprod(XW)
        XWX <- XWX + diag(1e-08, ncol(XWX))
        cholQ <- chol(XWX)
        M <- backsolve(cholQ, forwardsolve(t(cholQ), crossprod(Xj, wj * e)))
        M <- drop(M)

        ## Sample new parameters.
        b1 <- rmvnorm_cholQ(M, cholQ)

        ## Log-priors.
        p2 <- priors$p(b1)
        qbetaprop <- dmvnorm_cholQ(b1, M, cholQ)

        ## New fitted values.        
        fj <- drop(Xj %*% b1)

        ## Set up new predictor.
        eta[[j]] <- eta[[j]] + fj

        ## New parameters.
        peta <- family$map2par(eta)

        ## Compute new log likelihood.
        pibetaprop <- family$log_likelihood(par = peta, y = y)

        ## Compute new score and hessian.
        score <- deriv_checks(family$score[[j]](par = peta, y = y, id = j), is.weight = FALSE)
        hessian <- deriv_checks(family$hessian[[j]](par = peta, y = y, id = j), is.weight = TRUE)

        ## Weights.
        wj <- if(is.null(weights)) hessian else hessian * weights

        ## New working observations.
        z <- eta[[j]] + 1 / hessian * score

        ## New residuals.
        e <- z - eta2

        ## Compute mean and precision.
        XW <- Xj * sqrt(wj)
        XWX <- crossprod(XW)
        XWX <- XWX + diag(1e-08, ncol(XWX))
        cholQ <- chol(XWX)
        M <- backsolve(cholQ, forwardsolve(t(cholQ), crossprod(Xj, wj * e)))
        M <- drop(M)

        ## Log-priors.
        qbeta <- dmvnorm_cholQ(b0, M, cholQ)

        ## Acceptance probablity.
        alpha <- (pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1)

        ## Accept or reject?
        if(runif(1L) <= min(1, exp(alpha))) {
          fit[[j]]$coefficients <- b1
          fit[[j]]$fitted.values <- fj
        } else {
          eta[[j]] <- eta[[j]] - fj + fit[[j]]$fitted.values
        }

        ## Save.
        if(do_save) {
          samples[[j]]$p[isave, ] <- c(fit[[j]]$coefficients, min(1, exp(alpha)))
        }
      }

      ## Sample specials part.
      if(length(sterms[[j]])) {
        for(k in sterms[[j]]) {
          prop <- propose(x = specials[[k]], y = y, family = family,
            eta = eta, fitted = sfit[[j]][[k]], parameter = j,
            weights = weights, control = control)

          ## The smoothing variances are updated conditionally on the current
          ## coefficients before the Metropolis step. Keep these Gibbs/slice
          ## updates regardless of whether the coefficient proposal is
          ## accepted.
          sfit[[j]][[k]]$tau <- prop$tau
          sfit[[j]][[k]]$edf <- prop$edf
          sfit[[j]][[k]]$alpha <- prop$alpha

          ## Set new coefficient state.
          if(runif(1L) <= prop$alpha) {
            eta[[j]] <- eta[[j]] - sfit[[j]][[k]]$fitted.values + prop$fitted.values
            sfit[[j]][[k]]$fitted.values <- prop$fitted.values
            sfit[[j]][[k]]$coefficients <- prop$coefficients
          }

          ## Save.
          if(do_save) {
            samples[[j]]$s[[k]][isave, ] <- c(
              sfit[[j]][[k]]$coefficients,
              sfit[[j]][[k]]$tau, 
              sfit[[j]][[k]]$edf,
              if(is.null(sfit[[j]][[k]]$alpha)) 0 else sfit[[j]][[k]]$alpha
            )
          }
        }
      }

      ## Save eta.
      if(do_save) {
        track$eta[[j]] <- track$eta[[j]] + eta[[j]]
      }
    }

    ## Save global logLik / deviance once per saved iteration.
    if(do_save) {
      ll_iter <- family$log_likelihood(par = family$map2par(eta), y = y)
      track$logLik[isave] <- ll_iter
      track$deviance[isave] <- -2 * ll_iter
      isave <- isave + 1L
    }

    if(control$trace) {
      barfun(ptm, n.iter, iter, step, nstep)
    }
  }

  if(control$trace && interactive())
    cat("\n")

  ## Get mean coefficients.
  coef_lin <- list()
  for(j in np) {
    track$eta[[j]] <- track$eta[[j]] / nsave

    if(!is.null(samples[[j]]$p)) {
      coef_lin[[j]] <- apply(samples[[j]]$p, 2, mean, na.rm = TRUE)
      colnames(samples[[j]]$p) <- paste0(j, ".p.", colnames(samples[[j]]$p))
    }

    if(!is.null(samples[[j]]$s)) {
      for(k in names(samples[[j]]$s)) {
        nc <- ncol(specials[[k]]$X)

        colnames(samples[[j]]$s[[k]]) <- c(
          paste0(j, ".s.", k, ".", 1:nc),
          paste0(j, ".s.", k, ".tau", 1:length(specials[[k]]$S)),
          paste0(j, ".s.", k, ".edf"),
          paste0(j, ".s.", k, ".alpha")
        )

        kfit <- apply(samples[[j]]$s[[k]][, 1:nc, drop = FALSE], 1, function(b) {
          specials[[k]]$X %*% b
        })
        kfit <- apply(kfit, 1, mean)
        sfit[[j]][[k]]$fitted.values <- drop(kfit)
        lj <- grep(".tau", colnames(samples[[j]]$s[[k]]), fixed = TRUE)
        tau_samples <- samples[[j]]$s[[k]][, lj, drop = FALSE]
        sfit[[j]][[k]]$tau <- apply(tau_samples, 2, mean)
        sfit[[j]][[k]]$lambdas <- apply(1 / tau_samples, 2, mean)
        names(sfit[[j]][[k]]$lambdas) <- sub(
          ".tau", ".lambda", names(sfit[[j]][[k]]$lambdas), fixed = TRUE)
        lj <- grep(".edf", colnames(samples[[j]]$s[[k]]))
        sfit[[j]][[k]]$edf <- mean(samples[[j]]$s[[k]][, lj])
        lj <- grep(".alpha", colnames(samples[[j]]$s[[k]]))
        sfit[[j]][[k]]$alpha <- mean(samples[[j]]$s[[k]][, lj])
        sfit[[j]][[k]]$vcov <- cov(samples[[j]]$s[[k]][, 1:nc, drop = FALSE])
      }

      samples[[j]] <- cbind(samples[[j]]$p, do.call("cbind", samples[[j]]$s))
    } else {
      samples[[j]] <- do.call("cbind", samples[[j]])
    }
  }

  samples <- do.call("cbind", samples)

  ll <- family$log_likelihood(par = family$map2par(track$eta), y = y)

  Dbar <- mean(track$deviance, na.rm = TRUE)
  Dhat <- -2 * ll
  pD <- Dbar - Dhat
  DIC <- Dhat + 2 * pD

  dic <- list("Dbar" = Dbar, "Dhat" = Dhat, "pD" = pD, "DIC" = DIC)

  rval <- list(
    "fitted.values" = as.data.frame(track$eta),
    "fitted.specials" = sfit,
    "fitted.linear" = fit,
    "coefficients" = coef_lin,
    "iterations" = iter,
    "logLik" = ll, "control" = control,
    "nobs" = n,
    "deviance" = -2 * ll,
    "null.deviance" = dev0,
    "dev.reduction" = abs((dev0 - (-2 * ll)) / dev0),
    "dic" = dic,
    "nullmodel" = control$nullmodel,
    "samples" = samples
  )

  class(rval) <- c("bamlss2", "gamlss2")

  return(rval)
}

## Print info during sampling.
barfun <- function(ptm, n.iter, i, step, nstep, start = TRUE)
{
  ia <- interactive()
  if(i == 10 & start) {
    cat(if(ia) "\r" else "\n")
    elapsed <- c(proc.time() - ptm)[3]
    rt <- elapsed / i * (n.iter - i)
    rt <- if(rt > 60) {
      paste(formatC(format(round(rt / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(rt, 2), nsmall = 2), width = 5), "sec", sep = "")
    cat("|", rep(" ", nstep), "|   0% ", rt, sep = "")
    if(.Platform$OS.type != "unix" & ia) flush.console()
  }
  istep <- i %% step
  if(is.na(istep))
    istep <- 0
  if(istep == 0) {
    cat(if(ia) "\r" else "\n")
    p <- i / n.iter
    p <- paste("|", paste(rep("*", round(nstep * p)), collapse = ""),
      paste(rep(" ", round(nstep * (1 - p))), collapse = ""), "| ",
      formatC(round(p, 2) * 100, width = 3), "%", sep = "")
    elapsed <- c(proc.time() - ptm)[3]
    rt <- elapsed / i * (n.iter - i)
    rt <- if(rt > 60) {
      paste(formatC(format(round(rt / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(rt, 2), nsmall = 2), width = 5), "sec", sep = "")
    elapsed <- if(elapsed > 60) {
      paste(formatC(format(round(elapsed / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(elapsed, 2), nsmall = 2), width = 5), "sec", sep = "")
    cat(p, rt, elapsed, sep = " ")
    if(.Platform$OS.type != "unix" & ia) flush.console()
  }
}

## Generic log-prior function.
prior <- function(x, ...) {
  UseMethod("prior")
}

## Log-prior for mgcv smooth terms.
prior.mgcv.smooth <- function(x, ...)
{
  penalties <- x$S
  m <- length(penalties)

  a <- b <- 0.0001
  igs <- log((b^a)) - log(gamma(a))
  var_prior_fun <- function(tau) {
    sum(igs + (-a - 1) * log(tau) - b/tau)
  }

  ## A single scaled penalty has a fixed eigensystem:
  ## log|S/tau|+ = log|S|+ - rank(S) * log(tau).
  if(m == 1L) {
    ev <- eigen(penalties[[1L]], symmetric = TRUE, only.values = TRUE)$values
    tol <- max(ev) * 1e-12
    ev_pos <- ev[ev > tol]
    logdetS <- sum(log(ev_pos))
    rankS <- length(ev_pos)
  } else {
    ## The common null space is fixed for positive penalty weights.
    ## Restricting to its orthogonal complement turns every subsequent
    ## pseudo-determinant into an ordinary determinant.
    Ssum <- Reduce(`+`, penalties)
    ev <- eigen(Ssum, symmetric = TRUE)
    tol <- max(ev$values) * 1e-12
    U <- ev$vectors[, ev$values > tol, drop = FALSE]
    reduced_penalties <- lapply(penalties, function(S) {
      crossprod(U, S %*% U)
    })

    rank_aware_logdet <- function(tau) {
      P <- penalties[[1L]] / tau[1L]
      for(j in 2:m)
        P <- P + penalties[[j]] / tau[j]
      ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
      tol <- max(ev) * 1e-12
      ev_pos <- ev[ev > tol]
      sum(log(ev_pos))
    }

    ## Two penalties admit a one-time Demmler-Reinsch decomposition:
    ## every later log-determinant is then only a vector calculation.
    two_penalty_dr <- FALSE
    if(m == 2L && ncol(U)) {
      Rsum <- try(chol(reduced_penalties[[1L]] +
        reduced_penalties[[2L]]), silent = TRUE)
      if(!inherits(Rsum, "try-error")) {
        Ri <- backsolve(Rsum, diag(nrow(Rsum)))
        S1 <- crossprod(
          Ri,
          reduced_penalties[[1L]] %*% Ri
        )
        S1 <- (S1 + t(S1)) / 2
        penalty_eigen <- eigen(
          S1, symmetric = TRUE, only.values = TRUE
        )$values
        eig_tol <- sqrt(.Machine$double.eps)
        if(all(penalty_eigen > -eig_tol &
            penalty_eigen < 1 + eig_tol)) {
          penalty_eigen <- pmin(1, pmax(0, penalty_eigen))
          logdetSsum <- 2 * sum(log(diag(Rsum)))
          two_penalty_dr <- TRUE
        }
      }
    }
  }

  ## Slice updates repeatedly reuse both tau and the coefficient vector.
  cache <- new.env(parent = emptyenv())
  cache$tau <- NULL
  cache$logdetP <- NULL
  cache$gamma <- NULL
  cache$quadratics <- NULL

  function(parameters) {
    np <- length(parameters)
    i <- seq.int(np - m + 1L, np)
    nms <- names(parameters)

    ## The sampler stores tau at the end. Preserve support for externally
    ## supplied parameter vectors with a different named ordering.
    if(!is.null(nms) && !all(grepl(".tau", nms[i], fixed = TRUE))) {
      short_names <- vapply(
        strsplit(nms, ".s.", fixed = TRUE),
        function(z) z[length(z)],
        character(1L)
      )
      named_i <- grep(".tau", short_names, fixed = TRUE)
      if(length(named_i))
        i <- named_i
    }

    tau <- parameters[i]
    gamma_coef <- parameters[-i]
    ld <- var_prior_fun(tau)

    if(m == 1L) {
      Pgamma <- drop(penalties[[1L]] %*% gamma_coef)
      quad <- sum(gamma_coef * Pgamma) / tau[1L]
      logdetP <- logdetS - rankS * log(tau[1L])
    } else {
      tau_numeric <- as.numeric(tau)
      if(identical(tau_numeric, cache$tau)) {
        logdetP <- cache$logdetP
      } else {
        if(!ncol(U)) {
          logdetP <- 0
        } else if(two_penalty_dr) {
          penalty_values <- penalty_eigen / tau[1L] +
            (1 - penalty_eigen) / tau[2L]
          use_dr <- all(is.finite(penalty_values)) &&
            all(penalty_values > 0) &&
            max(tau) <= 1e08 * min(tau) &&
            min(penalty_values) > 1e-12 * max(penalty_values)
          logdetP <- if(use_dr) {
            logdetSsum + sum(log(penalty_values))
          } else {
            rank_aware_logdet(tau)
          }
        } else {
          Pr <- reduced_penalties[[1L]] / tau[1L]
          for(j in 2:m)
            Pr <- Pr + reduced_penalties[[j]] / tau[j]

          cholPr <- try(chol(Pr), silent = TRUE)
          use_chol <- !inherits(cholPr, "try-error")
          if(use_chol) {
            chol_diag <- abs(diag(cholPr))
            use_chol <- min(chol_diag) > 1e-06 * max(chol_diag)
          }

          logdetP <- if(use_chol) {
            2 * sum(log(chol_diag))
          } else {
            rank_aware_logdet(tau)
          }
        }

        cache$tau <- tau_numeric
        cache$logdetP <- logdetP
      }

      gamma_numeric <- as.numeric(gamma_coef)
      if(identical(gamma_numeric, cache$gamma)) {
        quadratics <- cache$quadratics
      } else {
        quadratics <- vapply(penalties, function(S) {
          sum(gamma_coef * drop(S %*% gamma_coef))
        }, numeric(1L))
        cache$gamma <- gamma_numeric
        cache$quadratics <- quadratics
      }
      quad <- sum(quadratics / tau)
    }

    lp <- 0.5 * logdetP - 0.5 * quad + ld
    lp[1L]
  }
}

## Generic propose function.
propose <- function(x, y, family, eta, fitted, parameter, weights = NULL, control = NULL)
{
  UseMethod("propose")
}

## Propose for mgcv smooth terms.
propose.mgcv.smooth <- function(x, y, family, eta, fitted,
  parameter, weights = NULL, control = NULL)
{
  ## Build chol(Q) and its mean. EDF is needed for the forward proposal only.
  build_QM_edf <- function(wj, e, tau, compute_edf = FALSE) {
    if(isTRUE(control$binning) && !is.null(x$binning)) {
      rw <- numeric(length(x$binning$nodups))
      rz <- numeric(length(x$binning$nodups))

      ## Reduce weights and response to unique rows.
      calc_Xe(x$binning$sorted.index, wj, e, rw, rz, x$binning$order)

      ## X'WX and X'We using reduced weights/response.
      XWX <- calc_XWX(x$X, 1/rw, x$sparse_index)
      XWz <- crossprod(x$X, rz)
    } else {
      ## The fused native path avoids materializing two weighted matrices.
      use_native <- ncol(x$X) >= 16L &&
        all(is.finite(wj)) && all(wj >= 0)
      if(use_native) {
        crossproducts <- calc_XWXz(x$X, wj, e)
        XWX <- crossproducts$XWX
        XWz <- crossproducts$XWz
      } else {
        XWX <- crossprod(x$X * sqrt(wj))
        XWz <- crossprod(x$X, wj * e)
      }
    }

    ## Keep the unpenalized crossproduct only for the forward EDF.
    if(compute_edf)
      XWX_unpenalized <- XWX

    ## Add penalties.
    for(jj in seq_along(tau))
      XWX <- XWX + 1/tau[jj] * x$S[[jj]]

    ## Stabilize and factorize.
    diag(XWX) <- diag(XWX) + 1e-08
    cholQ <- chol(XWX)

    ## Mean: M = Q^{-1}X'We, avoiding an explicit transpose.
    M <- backsolve(
      cholQ,
      backsolve(cholQ, XWz, transpose = TRUE)
    )
    M <- drop(M)

    ## tr(XW Q^{-1} XW') = tr(X'WX Q^{-1}). This avoids solving
    ## against all n rows of XW and is not needed for the reverse density.
    edf <- if(compute_edf) {
      sum(XWX_unpenalized * chol2inv(cholQ))
    } else {
      NULL
    }

    list("cholQ" = cholQ, "M" = M, "edf" = edf)
  }

  ## With parameter-independent Hessians, the reverse proposal has the
  ## same precision. Only its mean must then be updated.
  build_M <- function(wj, e, cholQ) {
    if(isTRUE(control$binning) && !is.null(x$binning)) {
      rw <- numeric(length(x$binning$nodups))
      rz <- numeric(length(x$binning$nodups))
      calc_Xe(x$binning$sorted.index, wj, e, rw, rz, x$binning$order)
      XWz <- crossprod(x$X, rz)
    } else {
      XWz <- crossprod(x$X, wj * e)
    }

    drop(backsolve(
      cholQ,
      backsolve(cholQ, XWz, transpose = TRUE)
    ))
  }

  ## Get parameters.
  peta <- family$map2par(eta)

  ## Compute old log-likelihood.
  pibeta <- family$log_likelihood(par = peta, y = y)

  ## Old parameters.
  b0 <- fitted$coefficients
  tau <- fitted$tau

  ## New shrinkage variance(s).
  if(!isTRUE(x$fixed)) {
    if(length(tau) > 1L) {
      theta <- c(b0, tau)
      tau_idx <- grep(".tau", names(theta), fixed = TRUE)
      for(jj in tau_idx) {
        theta <- uni.slice(theta, x, family,
          response = NULL, eta = NULL,
          id = parameter, j = jj,
          logPost = log_posterior,
          lower = 1e-08,
          log_likelihood = pibeta)
      }
      tau <- theta[tau_idx]
    } else {
      a <- x$rank / 2 + 0.0001
      b <- 0.5 * crossprod(b0, x$S[[1]]) %*% b0 + 0.0001
      nt <- names(tau)
      tau <- 1 / rgamma(1, a, b)
      names(tau) <- nt
    }
  }

  ## The variance update is a separate Gibbs/slice step. The coefficient
  ## Metropolis ratio is conditional on this updated variance, so both prior
  ## densities must be evaluated at the same value of tau.
  p1 <- x$prior(c(b0, tau))

  ## Working response and weights.
  ew <- .update(
    par = peta, y = y, eta = eta[[parameter]],
    family = family, which = parameter
  )
  z <- ew$eta
  hessian <- ew$weights

  ## Compute residuals.
  eta2 <- eta[[parameter]] <- eta[[parameter]] - fitted$fitted.values
  e <- z - eta2

  ## Weights.
  wj <- if(is.null(weights)) hessian else hessian * weights
  wj_forward <- wj

  ## Build the forward proposal, including EDF for saved output.
  tmp <- build_QM_edf(wj, e, tau, compute_edf = TRUE)
  cholQ <- tmp$cholQ
  M <- tmp$M
  edf <- tmp$edf
  cholQ_forward <- cholQ

  ## Sample new parameters and reuse the standardized draw for its density.
  draw <- rmvnorm_cholQ(M, cholQ, log_density = TRUE)
  b1 <- draw$sample
  qbetaprop <- draw$log_density

  ## Log-prior.
  p2 <- x$prior(c(b1, tau))

  ## New fitted values.
  fj0 <- drop(x$X %*% b1)
  fj <- if(isTRUE(control$binning) && !is.null(x$binning)) fj0[x$binning$match.index] else fj0

  ## Set up new predictor.
  eta[[parameter]] <- eta[[parameter]] + fj

  ## New parameters and log-likelihood.
  peta <- family$map2par(eta)
  pibetaprop <- family$log_likelihood(par = peta, y = y)

  ## New working response and weights.
  ew <- .update(
    par = peta, y = y, eta = eta[[parameter]],
    family = family, which = parameter
  )
  hessian <- ew$weights
  wj <- if(is.null(weights)) hessian else hessian * weights
  e <- ew$eta - eta2

  ## Reverse density: reuse Q when the Hessian weights are unchanged.
  if(identical(wj, wj_forward)) {
    cholQ <- cholQ_forward
    M <- build_M(wj, e, cholQ)
  } else {
    tmp <- build_QM_edf(wj, e, tau)
    cholQ <- tmp$cholQ
    M <- tmp$M
  }
  qbeta <- dmvnorm_cholQ(b0, M, cholQ)

  ## Acceptance probability.
  alpha <- (pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1)

  ## Assign names.
  names(b1) <- names(b0)

  ## Return state.
  fitted$fitted.values <- fj
  fitted$coefficients <- b1
  fitted$tau <- tau
  fitted$edf <- edf
  fitted$alpha <- min(1, exp(alpha))

  fitted
}

## Function to compute proportional log-posterior.
log_posterior <- function(coefficients, x, family, y,
  eta, parameter, log_likelihood = NULL)
{
  if(is.null(log_likelihood)) {
    eta[[parameter]] <- eta[[parameter]] + drop(x$X %*% coefficients[1:ncol(x$X)])
    log_likelihood <- family$log_likelihood(par = family$map2par(eta), y = y)
  }

  log_prior <- x$prior(coefficients)

  return(log_likelihood + log_prior)
}

rmvnorm_cholQ <- function(mean, cholQ, log_density = FALSE) {
  p <- length(mean)
  z <- rnorm(p)
  x <- mean + backsolve(cholQ, z, upper.tri = TRUE)

  if(!isTRUE(log_density))
    return(x)

  log_density <- sum(log(diag(cholQ))) -
    0.5 * sum(z * z) -
    0.5 * p * log(2 * pi)

  list("sample" = x, "log_density" = log_density)
}

dmvnorm_cholQ <- function(x, mean, cholQ) {
  p <- length(mean)
  r <- x - mean
  u <- cholQ %*% r
  quad <- sum(u * u)
  logdetQ <- 2 * sum(log(diag(cholQ)))
  0.5 * logdetQ - 0.5 * quad - 0.5 * p * log(2 * pi)
}

edf_from_cholQ_XP <- function(XW, cholQ) {
  B <- backsolve(cholQ, forwardsolve(t(cholQ), t(XW)))
  sum(XW * t(B))
}

## Univariate slice sampler.
uni.slice <- function(g, x, family, response, eta, id, j, ...,
  w = 1, m = 30, lower = -Inf, upper = +Inf, logPost)
{
  x0 <- g[j]
  gx0 <- logPost(g, x, family, response, eta, id, ...)

  ## Determine slice level (log).
  logy <- gx0 - rexp(1)

  ## Initial interval [L, R] of width w.
  u <- runif(1, 0, w)
  L <- x0 - u
  R <- x0 + (w - u)

  ## Step out.
  eval_at <- function(val) {
    old <- g[j]
    g[j] <- val
    out <- logPost(g, x, family, response, eta, id, ...)
    g[j] <- old
    out
  }

  if(is.infinite(m)) {
    repeat {
      if(L <= lower) break
      if(eval_at(L) <= logy) break
      L <- L - w
    }
    repeat {
      if(R >= upper) break
      if(eval_at(R) <= logy) break
      R <- R + w
    }
  } else if(m > 1) {
    J <- floor(runif(1, 0, m))
    K <- (m - 1) - J
    while(J > 0) {
      if(L <= lower) break
      if(eval_at(L) <= logy) break
      L <- L - w
      J <- J - 1
    }
    while(K > 0) {
      if(R >= upper) break
      if(eval_at(R) <= logy) break
      R <- R + w
      K <- K - 1
    }
  }

  ## Clamp to bounds.
  if(L < lower) L <- lower
  if(R > upper) R <- upper

  ## Shrinkage sampling.
  repeat {
    x1 <- runif(1, L, R)
    gx1 <- eval_at(x1)
    if(gx1 >= logy) {
      g[j] <- x1
      break
    }
    if(x1 > x0) R <- x1 else L <- x1
  }

  g
}

## Internal MCMC sampling function.
.mcmc <- function(x, y, specials, family, offsets, weights,
  start, xterms, sterms, control)
{
  if(is.null(control$trace))
    control$trace <- TRUE

  if(control$maxit[1] > 0) {
    if(control$trace[1L])
      cat(".. backfitting step\n")
    m <- RS(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
    start <- coef(m, full = TRUE, lambdas = TRUE)
  } else {
    stop("argument maxit must be > 1 for finding appropriate starting values!")
  }

  if(control$trace[1L])
    cat(".. MCMC step\n")

  m <- BS(x, y, specials, family, offsets, weights, start, xterms, sterms, control)

  return(m)
}

bamlss2 <- function(formula, n.iter = 1200, burnin = 200, thin = 1, maxit = 2, ...)
{
  call <- match.call()
  m <- call
  m[[1L]] <- as.name("gamlss2")
  m[["n.iter"]] <- n.iter
  m[["burnin"]] <- burnin
  m[["thin"]] <- thin
  m[["maxit"]] <- maxit
  m[["optimizer"]] <- getFromNamespace(".mcmc", "gamlss2")
  model <- eval(m, parent.frame())
  model$df <- model$dic$pD
  model$results <- results(model, data = model.frame(model))
  model$call <- call
  return(model)
}

mcmc <- function(object, n.iter = 1200, burnin = 200, thin = 1)
{
  if(!inherits(object, "gamlss2") && !inherits(object, "bamlss2")) {
    stop("wrong object supplied!")
  }

  ## Update control.
  object$control$n.iter <- n.iter
  object$control$burnin <- burnin
  object$control$thin <- thin

  ## Starting values (incl. lambdas).
  object$start <- coef(object, full = TRUE, lambdas = TRUE)

  ## Keep old samples (if any).
  samples0 <- object$samples

  ## Ensure x/y exist (works if model stored; like gamlss()).
  ## NOTE: if the object was fitted with control$light = TRUE and no model stored,
  ## you cannot reconstruct specials$X for MCMC. In that case, stop with message.
  if(is.null(object$y) || is.null(object$x)) {
    mf <- model.frame(object, keepresponse = TRUE)

    if(is.null(object$y)) {
      object$y <- model.response(mf)
      if(is.null(object$y)) {
        rn <- response_name(object$formula)
        object$y <- mf[, rn]
      }
    }

    if(is.null(object$x)) {
      object$x <- model.matrix(object, data = mf)
    }

    if(is.null(object$weights)) {
      object$weights <- model.weights(mf)
      if(!is.null(object$weights)) {
        if(length(object$weights) == 1L)
          object$weights <- rep.int(object$weights, nrow(mf))
        object$weights <- as.vector(object$weights)
        names(object$weights) <- rownames(mf)
      }
    }

    if(is.null(object$offsets)) {
      object$offsets <- model.offset(mf)
    }
  }

  ## If specials design matrices were dropped, BS can't run.
  if(!is.null(object$specials)) {
    for(k in seq_along(object$specials)) {
      if(is.null(object$specials[[k]][["X"]])) {
        stop("MCMC needs specials[[k]]$X, but it is NULL (probably fitted with control$light = TRUE). Refit with control$light = FALSE (and control$x=TRUE, control$y=TRUE).")
      }
    }
  }

  ## Build BS arguments explicitly (never drop weights/offsets).
  args <- list(
    x        = object$x,
    y        = object$y,
    specials = object$specials,
    family   = object$family,
    offsets  = object$offsets,
    weights  = object$weights,
    start    = object$start,
    xterms   = object$xterms,
    sterms   = object$sterms,
    control  = object$control
  )

  tstart <- proc.time()
  bs <- do.call(BS, args)
  elapsed <- as.numeric((proc.time() - tstart)["elapsed"])

  ## Merge BS output back into original object so terms/call/formula stay intact.
  for(nm in names(bs)) {
    object[[nm]] <- bs[[nm]]
  }

  ## Combine samples.
  if(!is.null(samples0) && nrow(samples0) > 0L) {
    object$samples <- rbind(samples0, object$samples)
  }

  ## Update derived summaries (these rely on terms being present).
  object$results <- results(object)
  object$df <- object$dic$pD
  object$elapsed <- elapsed
  object$call <- match.call()

  class(object) <- unique(c("bamlss2", class(object)))

  return(object)
}

## Propose and prior functions.
propose.elm <- propose.mgcv.smooth
prior.elm <- prior.mgcv.smooth

## Testing.
if(FALSE) {
  set.seed(123)

  n <- 1000

  d <- data.frame("x" = seq(-pi, pi, length = n))
  d$y <- 1.2 + sin(d$x) + rnorm(n, sd = exp(-1 + cos(d$x)))

  a <- gamlss2(y ~ s(x) | s(x), data = d)
  b <- mcmc(a)

  p <- predict(b)

  fit <- NULL
  for(j in c(0.025, 0.5, 0.975))
    fit <- cbind(fit, family(b)$quantile(p, j))

  par(mfrow = c(1, 2))
  plot(d)
  matplot(d$x, fit, type = "l", lty = 1, col = 4, lwd = 2, add = TRUE)
  matplot(b$samples, type = "l", lty = 1)
}

