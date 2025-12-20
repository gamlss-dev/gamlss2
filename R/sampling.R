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
  thinning <- control$thinning
  if(is.null(thinning))
    thinning <- 1L

  ## Type conversion.
  n.iter <- as.integer(n.iter)
  burnin <- as.integer(burnin)
  thinning <- as.integer(thinning)

  ## Numbers of samples to save.
  iterthin <- seq.int(burnin, n.iter, by = thinning)
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
      samples[[j]]$p <- matrix(NA, nrow = nsave, ncol = length(xterms[[j]]))
      colnames(samples[[j]]$p) <- xterms[[j]]
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
          if(!inherits(specials[[i]], "mgcv.smooth")) {
            stop("only mgcv smooth terms are allowed!")
          }
          sfit[[j]][[i]] <- list(
            "fitted.values" = rep(0.0, n),
            "edf" = 0.0,
            "coefficients" = setNames(rep(0.0, ncol(specials[[i]]$X)),
              colnames(specials[[i]]$X)),
            "tau" = setNames(rep(0.001, length(specials[[i]]$S)),
              paste0(".tau", seq_along(specials[[i]]$S)))
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
  dev0 <- -2 * family$logLik(y, family$map2par(eta))

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
        if(runif(1L) <= exp(alpha)) {
          fit[[j]]$coefficients <- b1
          fit[[j]]$fitted.values <- fj
        } else {
          eta[[j]] <- eta[[j]] - fj + fit[[j]]$fitted.values
        }

        ## Save.
        if(do_save) {
          samples[[j]]$p[isave, ] <- fit[[j]]$coefficients
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
    }

    if(do_save) {
      isave <- isave + 1L
    }

    if(control$trace) {
      barfun(ptm, n.iter, iter, step, nstep)
    }
  }

  if(control$trace && interactive())
    cat("\n")

  ## Get mean coefficients.
  coef_lin <- eta <- list()
  for(j in np) {
    eta[[j]] <- rep(0, n)

    if(!is.null(samples[[j]]$p)) {
      coef_lin[[j]] <- apply(samples[[j]]$p, 2, mean, na.rm = TRUE)
      colnames(samples[[j]]$p) <- paste0(j, ".p.", colnames(samples[[j]]$p))
      fit[[j]]$fitted.values <- drop(x[, xterms[[j]], drop = FALSE] %*% coef_lin[[j]])
      eta[[j]] <- eta[[j]] + fit[[j]]$fitted.values
    }

    if(!is.null(samples[[j]]$s)) {
      for(k in names(samples[[j]]$s)) {
        nc <- ncol(specials[[k]]$X)

        colnames(samples[[j]]$s[[k]]) <- c(
          paste0(j, ".s.", k, ".", 1:nc),
          paste0(j, ".s.", k, ".lambda", 1:length(specials[[k]]$S)),
          paste0(j, ".s.", k, ".edf"),
          paste0(j, ".s.", k, ".alpha")
        )

        cm <- apply(samples[[j]]$s[[k]][, 1:nc, drop = FALSE], 2, mean)
        sfit[[j]][[k]]$fitted.values <- drop(specials[[k]]$X %*% cm)
        lj <- grep(".lambda", colnames(samples[[j]]$s[[k]]))
        sfit[[j]][[k]]$lambda <- apply(samples[[j]]$s[[k]][, lj, drop = FALSE], 2, mean)
        lj <- grep(".edf", colnames(samples[[j]]$s[[k]]))
        sfit[[j]][[k]]$edf <- mean(samples[[j]]$s[[k]][, lj])
        lj <- grep(".alpha", colnames(samples[[j]]$s[[k]]))
        sfit[[j]][[k]]$alpha <- mean(samples[[j]]$s[[k]][, lj])
        sfit[[j]][[k]]$vcov <- cov(samples[[j]]$s[[k]][, 1:nc, drop = FALSE])

        eta[[j]] <- eta[[j]] + sfit[[j]][[k]]$fitted.values
      }

      samples[[j]] <- cbind(samples[[j]]$p, do.call("cbind", samples[[j]]$s))
    } else {
      samples[[j]] <- do.call("cbind", samples[[j]])
    }
  }

  samples <- do.call("cbind", samples)

  ll <- family$logLik(y, family$map2par(eta))

  rval <- list(
    "fitted.values" = as.data.frame(eta),
    "fitted.specials" = sfit,
    "fitted.linear" = fit,
    "coefficients" = coef_lin,
    "iterations" = iter,
    "logLik" = ll, "control" = control,
    "nobs" = length(eta[[1L]]),
    "deviance" = -2 * ll,
    "null.deviance" = dev0,
    "dev.reduction" = abs((dev0 - (-2 * ll)) / dev0),
    "nullmodel" = control$nullmodel,
    "samples" = samples
  )

  class(rval) <- c("gamlss2.mcmc", "gamlss2")

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
    i <- grep(".tau", names(parameters), fixed = TRUE)

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
  if(!x$fixed) {
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

  ## Compute mean and precision.
  XW <- x$X * sqrt(wj)
  XWX <- crossprod(XW)
  for(jj in seq_along(tau)) {
    XWX <- XWX + 1/tau[jj] * x$S[[jj]]
  }
  XWX <- XWX + diag(1e-08, ncol(XWX))
  cholQ <- chol(XWX)
  M <- backsolve(cholQ, forwardsolve(t(cholQ), crossprod(x$X, wj * e)))
  M <- drop(M)

  ## Degrees of freedom.
  edf <- edf_from_cholQ_XP(x$X, cholQ)

  ## Sample new parameters.
  b1 <- rmvnorm_cholQ(M, cholQ)

  ## Log-priors.
  p2 <- x$prior(c(b1, tau))
  qbetaprop <- dmvnorm_cholQ(b1, M, cholQ)

  ## New fitted values.        
  fj <- drop(x$X %*% b1)

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

  ## Compute mean and precision.
  XW <- x$X * sqrt(wj)
  XWX <- crossprod(XW)
  for(jj in seq_along(tau)) {
    XWX <- XWX + 1/tau[jj] * x$S[[jj]]
  }
  XWX <- XWX + diag(1e-08, ncol(XWX))
  cholQ <- chol(XWX)
  M <- backsolve(cholQ, forwardsolve(t(cholQ), crossprod(x$X, wj * e)))
  M <- drop(M)

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

edf_from_cholQ_XP <- function(X, cholQ) {
  B <- backsolve(cholQ, forwardsolve(t(cholQ), t(X)))
  sum(X * t(B))
}

## Univariate slice sampler.
uni.slice <- function(g, x, family, response, eta, id, j, ...,
  w = 1, m = 30, lower = -Inf, upper = +Inf, logPost)
{
  x0 <- g[j]
  gL <- gR <- g

  gx0 <- logPost(g, x, family, response, eta, id, ...)

  ## Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)

  ## Find the initial interval to sample from.
  u <- runif(1, 0, w)
  gL[j] <- g[j] - u
  gR[j] <- g[j] + (w - u)  ## should guarantee that g[j] is in [L, R], even with roundoff

  ## Expand the interval until its ends are outside the slice, or until
  ## the limit on steps is reached.
  if(is.infinite(m)) {
    repeat {
      if(gL[j] <= lower) break
      if(logPost(gL, x, family, response, eta, id, ...) <= logy) break
      gL[j] <- gL[j] - w
    }
    repeat {
      if(gR[j] >= upper) break
      if(logPost(gR, x, family, response, eta, id, ...) <= logy) break
      gR[j] <- gR[j] + w
    }
  } else {
    if(m > 1) {
      J <- floor(runif(1, 0, m))
      K <- (m - 1) - J
      while(J > 0) {
        if(gL[j] <= lower) break
        if(logPost(gL, x, family, response, eta, id, ...) <= logy) break
        gL[j] <- gL[j] - w
        J <- J - 1
      }
      while(K > 0) {
        if(gR[j] >= upper) break
        if(logPost(gR, x, family, response, eta, id, ...) <= logy) break
        gR[j] <- gR[j] + w
        K <- K - 1
      }
    }
  }

  ## Shrink interval to lower and upper bounds.
  if(gL[j] < lower) {
    gL[j] <- lower
  }
  if(gR[j] > upper) {
    gR[j] <- upper
  }

  ## Sample from the interval, shrinking it on each rejection.
  repeat {
    g[j] <- runif(1, gL[j], gR[j])

    gx1 <- logPost(g, x, family, response, eta, id, ...)

    if(gx1 >= logy) break

    if(g[j] > x0) {
      gR[j] <- g[j]
    } else {
      gL[j] <- g[j]
    }
  }

  ## Return the point sampled
  return(g)
}

## Testing.
if(FALSE) {
  set.seed(123)

  n <- 1000

  d <- data.frame("x" = seq(-pi, pi, length = n))
  d$y <- 1.2 + sin(d$x) + rnorm(n, sd = exp(-1 + cos(d$x)))

  m <- gamlss2(y ~ s(x,k=40) | s(x,k=40), data = d, optimizer = RS, maxit = 1)

  cm <- coef(m, full = TRUE, lambdas = TRUE)

  b <- gamlss2(y ~ s(x,k=40) | s(x,k=40), data = d, optimizer = BS, start = cm)

  p <- predict(b, FUN = median)

  fit <- NULL
  for(j in c(0.025, 0.5, 0.975))
    fit <- cbind(fit, family(b)$quantile(j, p))

  par(mfrow = c(1, 2))
  plot(d)
  matplot(d$x, fit, type = "l", lty = 1, col = 4, lwd = 2, add = TRUE)
  matplot(b$samples, type = "l", lty = 1)
}

