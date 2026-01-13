## Function to compute multivariate normal samples
## given the maximum likelihood estimator.
sampling <- function(object, R = 100, ...)
{
  V <- vcov(object, ...)
  ok <- TRUE
  Cv <- tryCatch(chol(V), error=function(e) { ok <<- FALSE; NULL })
  if(!ok) {
    V <- as.matrix(Matrix::nearPD(V)$mat)
    Cv <- chol(V)
  }
  cb <- coef(object, dropall = FALSE, ...)
  sc <- rnorm(R * length(cb))
  sc <- t(cb + t(Cv) %*% matrix(sc, nrow = length(cb), ncol = R))
  d <- drop(cb - apply(sc, 2, mean))
  sc <- t(t(sc) + d)
  colnames(sc) <- names(cb)
  return(sc)
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
              paste0(j, ".s.", i, ".", seq_along(ncol(specials[[i]]$X)))),
            "tau" = setNames(rep(0.001, length(specials[[i]]$S)),
              paste0(j, ".s.", i, ".tau", seq_along(specials[[i]]$S)))
          )
          samples[[j]]$s[[i]] <- matrix(NA, nrow = nsave,
            ncol = ncol(specials[[i]]$X) + length(specials[[i]]$S) + 2L)
          if(!is.null(cstart)) {
            sj <- grep(paste0(j, ".s.", i), names(cstart), fixed = TRUE, value = TRUE)
            sjb <- sj[!grepl(".lambda", sj, fixed = TRUE)]
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

  ## Null deviance.
  dev0 <- -2 * family$logLik(y, family$map2par(etastart))

  ## Estimate intercept only model first.
  if(isTRUE(control$nullmodel) & length(xterms)) {
    beta <- ieta <- list()
    for(j in np) {
      beta[[j]] <- as.numeric(fit[[j]]$coefficients["(Intercept)"])
      ieta[[j]] <- rep(beta[[j]], n)
    }
    beta <- unlist(beta)

    if(!any(is.na(beta))) {
      lli <- family$logLik(y, family$map2par(ieta))

      fn_ll <- function(par) {
        for(j in np)
          ieta[[j]] <- rep(par[j], n)
        ll <- family$logLik(y, family$map2par(ieta)) - 1e-05 * sum(par^2)
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
        dev0 <- -2 * family$logLik(y, family$map2par(eta))
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
      if(control$fixed[[j]])
        stop("fixed parameters not supported yet!")

      ## Sampling linear part.
      if(length(xterms[[j]])) {
        ## Get parameters.
        peta <- family$map2par(eta)

        ## Compute old log-likelihood.
        pibeta <- family$logLik(y, peta)

        ## Old parameters.
        b0 <- fit[[j]]$coefficients

        ## Log-prior.
        p1 <- priors$p(b0)

        ## Derivatives.
        score <- deriv_checks(family$score[[j]](y, peta, id = j), is.weight = FALSE)
        hess <- deriv_checks(family$hess[[j]](y, peta, id = j), is.weight = TRUE)

        ## Working response.
        z <- eta[[j]] + 1 / hess * score

        ## Compute residuals.
        eta2 <- eta[[j]] <- eta[[j]] - fit[[j]]$fitted.values
        e <- z - eta2

        ## Weights.
        wj <- if(is.null(weights)) hess else hess * weights

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
        pibetaprop <- family$logLik(y, peta)

        ## Compute new score and hess.
        score <- deriv_checks(family$score[[j]](y, peta, id = j), is.weight = FALSE)
        hess <- deriv_checks(family$hess[[j]](y, peta, id = j), is.weight = TRUE)

        ## Weights.
        wj <- if(is.null(weights)) hess else hess * weights

        ## New working observations.
        z <- eta[[j]] + 1 / hess * score

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

          ## Set new state.
          if(runif(1L) <= prop$alpha) {
            eta[[j]] <- eta[[j]] - sfit[[j]][[k]]$fitted.values + prop$fitted.values
            sfit[[j]][[k]]$fitted.values <- prop$fitted.values
            sfit[[j]][[k]]$coefficients <- prop$coefficients
            sfit[[j]][[k]]$tau <- prop$tau
            sfit[[j]][[k]]$edf <- prop$edf
            sfit[[j]][[k]]$alpha <- prop$alpha
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
      ll_iter <- family$logLik(y, family$map2par(eta))
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
          paste0(j, ".s.", k, ".lambdas", 1:length(specials[[k]]$S)),
          paste0(j, ".s.", k, ".edf"),
          paste0(j, ".s.", k, ".alpha")
        )

        kfit <- apply(samples[[j]]$s[[k]][, 1:nc, drop = FALSE], 1, function(b) {
          specials[[k]]$X %*% b
        })
        kfit <- apply(kfit, 1, mean)
        sfit[[j]][[k]]$fitted.values <- drop(kfit)
        lj <- grep(".lambdas", colnames(samples[[j]]$s[[k]]))
        sfit[[j]][[k]]$lambdas <- apply(samples[[j]]$s[[k]][, lj, drop = FALSE], 2, mean)
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

  ll <- family$logLik(y, family$map2par(track$eta))

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
  function(parameters) {
    nms <- names(parameters)
    i <- integer(0)
    if(!is.null(nms)) {
      nms <- sapply(strsplit(nms, ".s."), function(x) x[length(x)])
      i <- grep(".tau", nms, fixed = TRUE)
    }
    if(length(i) == 0L) {
      m <- length(x$S)
      i <- (length(parameters) - m + 1L):length(parameters)
    }

    tau <- parameters[i]
    gamma <- parameters[-i]

    a <- b <- 0.0001
    igs <- log((b^a)) - log(gamma(a))
    var_prior_fun <- function(tau) {
      igs + (-a - 1) * log(tau) - b/tau
    }

    ld <- 0
    P <- 0

    for(j in seq_along(tau)) {
      P <- P + 1/tau[j] * x$S[[j]]
      ld <- ld + var_prior_fun(tau[j])
    }

    if(is.null(dim(P))) {
      P <- matrix(P, 1L, 1L)
    }

    ## Pseudo log-determinant (rank-aware).
    ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
    tol <- max(ev) * 1e-12
    ev_pos <- ev[ev > tol]
    logdetP <- sum(log(ev_pos))

    ## Quadratic form.
    quad <- drop(crossprod(gamma, P %*% gamma))

    lp <- 0.5 * logdetP - 0.5 * quad + ld

    return(lp[1L])
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
  ## Helper to build chol(Q), mean M, and edf using (optional) binning.
  build_QM_edf <- function(wj, e, tau) {
    if(isTRUE(control$binning) && !is.null(x$binning)) {
      rw <- numeric(length(x$binning$nodups))
      rz <- numeric(length(x$binning$nodups))

      ## Reduce weights and response to unique rows.
      calc_Xe(x$binning$sorted.index, wj, e, rw, rz, x$binning$order)

      ## X'WX and X'We using reduced weights/response.
      XWX <- calc_XWX(x$X, 1/rw, x$sparse_index)
      XWz <- crossprod(x$X, rz)

      ## For edf we need W^{1/2}X on unique rows.
      XW <- x$X * sqrt(rw)
    } else {
      ## Full data.
      XWX <- crossprod(x$X * sqrt(wj))
      XWz <- crossprod(x$X, wj * e)
      XW  <- x$X * sqrt(wj)
    }

    ## Add penalties.
    for(jj in seq_along(tau)) {
      XWX <- XWX + 1/tau[jj] * x$S[[jj]]
    }

    ## Stabilize and factorize.
    XWX <- XWX + diag(1e-08, ncol(XWX))
    cholQ <- chol(XWX)

    ## Mean: M = Q^{-1}X'We (no explicit inverse).
    M <- backsolve(cholQ, forwardsolve(t(cholQ), XWz))
    M <- drop(M)

    ## EDF: tr( XW Q^{-1} XW' ).
    edf <- edf_from_cholQ_XP(XW, cholQ)

    list("cholQ" = cholQ, "M" = M, "edf" = edf)
  }

  ## Get parameters.
  peta <- family$map2par(eta)

  ## Compute old log-likelihood.
  pibeta <- family$logLik(y, peta)

  ## Old parameters.
  b0 <- fitted$coefficients
  tau <- fitted$tau

  ## Log-prior.
  p1 <- x$prior(c(b0, tau))

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
      a <- x$rank / 2 + 0.001
      b <- 0.5 * crossprod(b0, x$S[[1]]) %*% b0 + 0.001
      nt <- names(tau)
      tau <- 1 / rgamma(1, a, b)
      names(tau) <- nt
    }
  }

  ## Derivatives.
  score <- deriv_checks(
    family$score[[parameter]](y, peta, id = parameter),
    is.weight = FALSE
  )
  hess <- deriv_checks(
    family$hess[[parameter]](y, peta, id = parameter),
    is.weight = TRUE
  )

  ## Working response.
  z <- eta[[parameter]] + 1 / hess * score

  ## Compute residuals.
  eta2 <- eta[[parameter]] <- eta[[parameter]] - fitted$fitted.values
  e <- z - eta2

  ## Weights.
  wj <- if(is.null(weights)) hess else hess * weights

  ## Build proposal precision + mean (+ edf) using binning-aware code.
  tmp <- build_QM_edf(wj, e, tau)
  cholQ <- tmp$cholQ
  M <- tmp$M
  edf <- tmp$edf

  ## Sample new parameters.
  b1 <- rmvnorm_cholQ(M, cholQ)

  ## Log-priors.
  p2 <- x$prior(c(b1, tau))
  qbetaprop <- dmvnorm_cholQ(b1, M, cholQ)

  ## New fitted values.
  fj0 <- drop(x$X %*% b1)
  fj <- if(isTRUE(control$binning) && !is.null(x$binning)) fj0[x$binning$match.index] else fj0

  ## Set up new predictor.
  eta[[parameter]] <- eta[[parameter]] + fj

  ## New parameters.
  peta <- family$map2par(eta)

  ## Compute new log likelihood.
  pibetaprop <- family$logLik(y, peta)

  ## Compute new score and hess.
  score <- deriv_checks(
    family$score[[parameter]](y, peta, id = parameter),
    is.weight = FALSE
  )
  hess <- deriv_checks(
    family$hess[[parameter]](y, peta, id = parameter),
    is.weight = TRUE
  )

  ## Weights.
  wj <- if(is.null(weights)) hess else hess * weights

  ## New working observations.
  z <- eta[[parameter]] + 1 / hess * score

  ## New residuals.
  e <- z - eta2

  ## Reverse density: rebuild Q,M at the new state.
  tmp <- build_QM_edf(wj, e, tau)
  cholQ <- tmp$cholQ
  M <- tmp$M

  ## Log-priors.
  qbeta <- dmvnorm_cholQ(b0, M, cholQ)

  ## Acceptance probablity.
  alpha <- (pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1)

  ## Assign names.
  names(b1) <- names(b0)

  ## Return state.
  fitted$fitted.values <- fj
  fitted$coefficients <- b1
  fitted$tau <- tau
  fitted$edf <- edf
  fitted$alpha <- min(1, exp(alpha))

  return(fitted)
}

## Function to compute proportional log-posterior.
log_posterior <- function(coefficients, x, family, y,
  eta, parameter, log_likelihood = NULL)
{
  if(is.null(log_likelihood)) {
    eta[[parameter]] <- eta[[parameter]] + drop(x$X %*% coefficients[1:ncol(x$X)])
    log_likelihood <- family$logLik(y, family$map2par(eta))
  }

  log_prior <- x$prior(coefficients)

  return(log_likelihood + log_prior)
}

rmvnorm_cholQ <- function(mean, cholQ) {
  p <- length(mean)
  z <- rnorm(p)
  x <- mean + backsolve(cholQ, z, upper.tri = TRUE)
  x
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
  object$df <- get_df(object)
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
    fit <- cbind(fit, family(b)$quantile(j, p))

  par(mfrow = c(1, 2))
  plot(d)
  matplot(d$x, fit, type = "l", lty = 1, col = 4, lwd = 2, add = TRUE)
  matplot(b$samples, type = "l", lty = 1)
}

