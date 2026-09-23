## Rigby and Stasinopoulos algorithm.
RS <- function(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
{
  ## Number of observations.
  n <- if(is.null(dim(y))) length(y) else nrow(y)

  ## Parameter names. FIXME: TRUE/FALSE?
  np <- family$names

  ## Caches belong to this fit only. Only unchanged, generated distribution
  ## callbacks opt in; user replacements may depend on external state.
  use.cache <- !identical(control$rs.cache, FALSE)
  cached.functions <- attr(family, "rs.cache", exact = TRUE)
  map2par <- family$map2par
  log_likelihood <- family$log_likelihood
  pdf <- family$pdf
  if(use.cache) {
    if(identical(map2par, cached.functions$map2par))
      map2par <- rs_cached_function(map2par)
    if(identical(log_likelihood, cached.functions$log_likelihood))
      log_likelihood <- rs_cached_function(log_likelihood)
    if(identical(pdf, cached.functions$pdf))
      pdf <- rs_cached_function(pdf)
  }
  linear.cache <- smooth.cache <- list()
  if(use.cache) {
    for(j in np) {
      linear.cache[[j]] <- new.env(parent = emptyenv())
      smooth.cache[[j]] <- list()
      for(k in sterms[[j]])
        smooth.cache[[j]][[k]] <- new.env(parent = emptyenv())
    }
  }

  ## Stepwise candidate fits only need the objective and degrees of freedom.
  ## Expensive inferential output is computed by
  ## the final, full fit.
  stepwise_candidate <- isTRUE(control$.stepwise_candidate)

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

  ## Trial updates can temporarily leave the parameter space, especially for
  ## positive parameters fitted with an identity link.  Such a trial has
  ## likelihood -Inf and must be rejected by the existing step-length logic;
  ## it is not a fatal family error.  Likelihood evaluations of the current
  ## accepted fit remain unprotected so genuine family errors are visible.
  candidate_log_likelihood <- function(eta) {
    tryCatch({
      par <- map2par(eta)
      if(is.null(weights)) {
        log_likelihood(par = par, y = y)
      } else {
        sum(pdf(par = par, y = y, log = TRUE) * weights, na.rm = TRUE)
      }
    }, error = function(e) -Inf)
  }

  ## Set control parameters.
  ## Stopping criterion.
  eps <- control$eps
  if(is.null(eps))
    eps <- 0.00001 ## sqrt(.Machine$double.eps)
  if(!is.numeric(eps) || !length(eps) || length(eps) > 3L ||
      any(!is.finite(eps)) || any(eps <= 0))
    stop("argument eps must contain one to three positive finite values!")
  if(length(eps) < 2L)
    eps <- c(eps, eps)
  if(length(eps) < 3L)
    eps <- c(eps, eps[2L])
  stop.eps <- eps
  control$eps <- stop.eps
  eps <- eps + 1

  ## The step length control.
  if(is.null(control$step))
    control$step <- 1
  if(!is.numeric(control$step) || length(control$step) != 1L ||
      !is.finite(control$step))
    stop("argument step must be one finite value!")
  if((control$step > 1) | (control$step < 0))
    control$step <- 1
  if(is.null(control$autostep))
    control$autostep <- TRUE
  else
    control$autostep <- isTRUE(control$autostep)

  if(is.null(control$sigma.tol))
    control$sigma.tol <- FALSE
  if(identical(control$sigma.tol, FALSE))
    control$sigma.tol <- 0
  if(!is.numeric(control$sigma.tol) || length(control$sigma.tol) != 1L ||
      !is.finite(control$sigma.tol) || control$sigma.tol < 0 ||
      control$sigma.tol >= 1)
    stop("argument sigma.tol must be one number between zero and one!")

  ## Maximum number of backfitting iterations.
  maxit <- control$maxit
  if(is.null(maxit))
    maxit <- 20L
  if(!is.numeric(maxit) || !length(maxit) || length(maxit) > 3L ||
      any(!is.finite(maxit)) || any(maxit < 1) ||
      any(maxit != floor(maxit)) || any(maxit > .Machine$integer.max))
    stop("argument maxit must contain one to three positive integers!")
  maxit <- as.integer(maxit)
  if(length(maxit) < 2L)
    maxit <- c(maxit, 50L)
  if(length(maxit) < 3L)
    maxit <- c(maxit, maxit[2L])
  control$maxit <- maxit

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
          fit[[j]]$fitted.values <- drop(x[, names(fit[[j]]$coefficients), drop = FALSE] %*% fit[[j]]$coefficients)
          nes[[j]] <- TRUE
        }
      }
      eta[[j]] <- fit[[j]]$fitted.values
    }
    if(length(sterms)) {
      if(length(sterms[[j]])) {
        sfit[[j]] <- list()
        for(i in sterms[[j]]) {
          sfit[[j]][[i]] <- list("fitted.values" = rep(0.0, n), "edf" = 0.0, "selected" = FALSE)
          if(!is.null(cstart)) {
            prefix <- paste0(j, ".s.", i, ".")
            sj <- names(cstart)[startsWith(names(cstart), prefix)]
            sjb <- sj[!grepl(".lambda", sj, fixed = TRUE)]
            sjb <- sjb[!grepl(".tau", sjb, fixed = TRUE)]
            if(length(sjb)) {
              if(!is.null(specials[[i]]$X) &&
                  length(sjb) == ncol(specials[[i]]$X)) {
                bstart <- as.numeric(cstart[sjb])
                sfit[[j]][[i]]$fitted.values <- drop(specials[[i]]$X %*% bstart)
                sfit[[j]][[i]]$coefficients <- bstart
                sfit[[j]][[i]]$vcov <- matrix(
                  NA_real_, length(bstart), length(bstart)
                )
                sfit[[j]][[i]]$transfer <- list(
                  "coefficients" = bstart,
                  "names" = colnames(specials[[i]]$X)
                )
                sjl <- sj[grepl(".lambda", sj, fixed = TRUE)]
                if(length(sjl))
                  sfit[[j]][[i]]$transfer$lambdas <- as.numeric(cstart[sjl])
                sfit[[j]][[i]]$.from_start <- TRUE
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

  if(!is.null(control$fixed)) {
    for(j in np) {
      if(control$fixed[[j]] && is.null(start[[j]])) {
        link <- make.link2(family$links[[j]])
        fit[[j]]$coefficients["(Intercept)"] <- link$linkfun(control$fixed[[j]])
        eta[[j]] <- rep(fit[[j]]$coefficients["(Intercept)"], n)
        fit[[j]]$fitted.values <- eta[[j]]
        etastart[[j]] <- eta[[j]]
      }
    }
  }

  ## Null deviance.
  dev0 <- -2 * log_likelihood(par = map2par(etastart), y = y)

  ## Estimate intercept only model first.
  run_nullmodel <- !stepwise_candidate ||
    (isTRUE(control$initialize) && is.null(start))
  if(run_nullmodel && isTRUE(control$nullmodel) & length(unlist(xterms))) {
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
      lli <- log_likelihood(par = map2par(ieta), y = y)

      fn_ll <- function(par) {
        for(j in np) {
          if(control$fixed[[j]])
            par[j] <- beta[j]
          ieta[[j]] <- rep(par[j], n)
          if(!is.null(offsets)) {
            if(!is.null(offsets[[j]]))
              ieta[[j]] <- ieta[[j]] + offsets[[j]]
          }
        }
        ll <- log_likelihood(par = map2par(ieta), y = y) - lambda * sum(par^2)
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
  if(length(family$hessian) < 2L)
    CGk <- Inf
  if(!any(grepl(":", names(family$hessian), fixed = TRUE)))
    CGk <- Inf
  if(is.finite(CGk))
    eta_old <- eta
  ## A global RS iteration is already one complete sweep over all parameters.
  ## Repeating that sweep here bypasses the global convergence check and keeps
  ## the outer iteration counter fixed.  Only CG needs internal sweeps because
  ## its working quantities are held fixed while solving the coupled system.
  maxit_RS <- 1L
  maxit_CG <- maxit[3L]

  ## Track iterations
  iter <- c(0L, 0L)
  dev.warn <- 0L
  dev.warn.max <- 0

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

  ## Interpolate a fitted term contribution.  For fixed-design terms automatic
  ## step length control must not leave the fitted values, coefficients, and
  ## the warm start at different points.
  blend_fit <- function(old, new, alpha, coefficients = TRUE) {
    rval <- new
    if(coefficients && !is.null(new$coefficients)) {
      b0 <- old$coefficients
      if(is.null(b0))
        b0 <- rep(0, length(new$coefficients))
      if(length(b0) != length(new$coefficients) ||
          !identical(dim(b0), dim(new$coefficients)))
        stop("cannot interpolate fitted term coefficients!")
      rval$coefficients <- alpha * new$coefficients + (1 - alpha) * b0
      if(!is.null(rval$transfer))
        rval$transfer$coefficients <- rval$coefficients
    }
    if(!is.null(new$fitted.values)) {
      f0 <- old$fitted.values
      if(is.null(f0))
        f0 <- rep(0, length(new$fitted.values))
      rval$fitted.values <- alpha * new$fitted.values + (1 - alpha) * f0
    }
    rval
  }

  blend_parameter <- function(old.fit, new.fit, old.sfit, new.sfit, alpha) {
    rval <- list("fit" = blend_fit(old.fit, new.fit, alpha), "sfit" = new.sfit)
    if(length(new.sfit)) {
      for(k in names(new.sfit)) {
        sk <- specials[[k]]
        fixed.design <- inherits(sk, c("mgcv.smooth", "X %*% b")) ||
          (!inherits(sk, c("smooth", "special")) &&
            is.null(sk$special.wfit))
        rval$sfit[[k]] <- blend_fit(old.sfit[[k]], new.sfit[[k]], alpha,
          fixed.design)
      }
    }
    rval
  }

  ## Smoothing parameters and coefficients form one state. If a criterion
  ## selects a new smoothing parameter, accept or reject that state as a
  ## whole; interpolating only its coefficients leaves an invalid warm start.
  smoothing_update <- function(sfit) {
    if(length(intersect(names(sfit), np)))
      sfit <- unlist(sfit, recursive = FALSE)
    state <- unlist(lapply(sfit, function(x) list(x$.rs_smoothing)),
      recursive = FALSE)
    state <- Filter(function(x) is.list(x) && isTRUE(x$changed), state)
    if(!length(state))
      return(list("changed" = FALSE, "scored" = FALSE,
        "improved" = FALSE))
    scored <- all(vapply(state, function(x) isTRUE(x$scored), logical(1L)))
    improved <- scored &&
      all(vapply(state, function(x) isTRUE(x$improved), logical(1L)))
    list("changed" = TRUE, "scored" = scored, "improved" = improved)
  }

  smooth_penalty <- function(old, new, alpha) {
    if(!length(new))
      return(0)
    if(length(intersect(names(new), np))) {
      return(sum(vapply(names(new), function(j)
        smooth_penalty(old[[j]], new[[j]], alpha), numeric(1L))))
    }

    value <- 0
    for(k in names(new)) {
      b1 <- new[[k]]$coefficients
      if(is.null(b1))
        next
      b0 <- old[[k]]$coefficients
      if(is.null(b0))
        b0 <- rep(0, length(b1))
      if(length(b0) != length(b1))
        return(Inf)
      b <- alpha * b1 + (1 - alpha) * b0
      P <- new[[k]]$penalty
      if(is.null(P)) {
        S <- specials[[k]][["S", exact = TRUE]]
        if(is.null(S) || !length(S) || is.null(new[[k]]$lambdas))
          next
        P <- matrix(0, length(b), length(b))
        lambda <- rep(new[[k]]$lambdas, length.out = length(S))
        for(i in seq_along(S))
          P <- P + lambda[i] * S[[i]]
      }
      if(!identical(dim(P), c(length(b), length(b))))
        next
      value <- value + drop(crossprod(b, P %*% b))
    }
    value
  }

  ## Backtrack a complete RS parameter update or a complete CG correction
  ## sweep.  A trial outside the parameter space has likelihood -Inf.
  find_step <- function(eta0, eta1, ll0, initial, old.sfit, new.sfit) {
    if(initial <= 0)
      return(list("eta" = eta0, "logLik" = ll0, "step" = 0, "accepted" = FALSE))
    alpha <- initial
    penalty0 <- smooth_penalty(old.sfit, new.sfit, 0)
    objective0 <- ll0 - 0.5 * penalty0
    tolerance <- sqrt(.Machine$double.eps) * (1 + abs(objective0))
    repeat {
      etai <- Map(function(a, b) a + alpha * (b - a), eta0, eta1)
      names(etai) <- names(eta0)
      ll1 <- candidate_log_likelihood(etai)
      penalty1 <- smooth_penalty(old.sfit, new.sfit, alpha)
      objective1 <- ll1 - 0.5 * penalty1
      if(is.finite(objective1) && objective1 >= objective0 - tolerance) {
        return(list("eta" = etai, "logLik" = ll1,
          "step" = alpha, "accepted" = TRUE,
          "objective" = penalty0 != 0 || penalty1 != 0))
      }
      if(!control$autostep || alpha <= sqrt(.Machine$double.eps))
        break
      alpha <- alpha * 0.5
    }
    list("eta" = eta0, "logLik" = ll0, "step" = 0, "accepted" = FALSE)
  }

  safeguard <- function(eta0, eta1, ll0, initial, old.sfit, new.sfit) {
    smoothing <- smoothing_update(new.sfit)
    if(!smoothing$changed)
      return(find_step(eta0, eta1, ll0, initial, old.sfit, new.sfit))

    ll1 <- candidate_log_likelihood(eta1)
    accept <- is.finite(ll1) && if(smoothing$scored)
      smoothing$improved else TRUE
    if(accept) {
      return(list("eta" = eta1, "logLik" = ll1, "step" = 1,
        "accepted" = TRUE, "objective" = TRUE))
    }
    list("eta" = eta0, "logLik" = ll0, "step" = 0,
      "accepted" = FALSE, "objective" = FALSE)
  }

  set_step <- function(j, alpha) {
    if(length(xterms[[j]]))
      step[[j]]$xterms <<- alpha
    if(length(sterms[[j]]))
      step[[j]]$sterms[] <<- alpha
  }

  set_all_steps <- function(alpha) {
    for(j in np)
      if(!control$fixed[[j]])
        set_step(j, alpha)
  }

  ## Start outer loop.
  while((eps[1L] > stop.eps[1L]) && (iter[1L] < maxit[1L])) {
    objective.accepted <- FALSE
    safeguard.failed <- FALSE
    ## Old log-likelihood.
    if(is.null(weights)) {
      llo0 <- log_likelihood(par = map2par(eta), y = y)
    } else {
      llo0 <- sum(pdf(par = map2par(eta), y = y, log = TRUE) * weights, na.rm = TRUE)
    }

    ## CG = k means that the kth displayed outer iteration uses CG.  The
    ## counter itself is zero-based at the top of this loop.
    use_CG <- iter[1L] + 1L >= CGk

    ## For CG.
    if(use_CG) {
      eta_old <- if(iter[1L] > 0L) eta else etastart
      par <- if(iter[1L] > 0L) {
        map2par(eta)
      } else {
        map2par(etastart)
      }
      ew_CG <- list()
      for(j in np) {
        ew_CG[[j]] <- .update(par = par, y = y,
          eta = if(iter[1L] > 0L) eta[[j]] else etastart[[j]],
          family = family, which = j)
      }
    }

    eps_outer <- 1
    iter_outer <- 0

    while((eps_outer > stop.eps[3L]) &&
        (iter_outer < if(use_CG) maxit_CG else maxit_RS)) {
      if(use_CG) {
        if(is.null(weights)) {
          outer_ll0 <- log_likelihood(par = map2par(eta), y = y)
        } else {
          outer_ll0 <- sum(pdf(par = map2par(eta), y = y, log = TRUE) * weights, na.rm = TRUE)
        }
        eta_sweep <- eta
        fit_sweep <- fit
        sfit_sweep <- sfit
        penalty_sweep <- penalty
      }

      for(j in np) {
        ## Check if paramater is fixed.
        if(control$fixed[[j]])
          next

        ## Outer loop working response and weights.
        par <- if(iter[1L] > 0L) {
          map2par(eta)
        } else {
          map2par(etastart)
        }

        ## Compute working response z and weights hessian from family.
        ## Cole and Green adjustment.
        if(use_CG) {
          h <- grep(paste0(j, ":"), names(family$hessian), value = TRUE)
          if(length(h)) {
            adj <- 0.0
            for(l in seq_along(h)) {
              parts <- strsplit(h[l], ":", fixed = TRUE)[[1]]
              k <- parts[2L]
              hessian_l <- family$hessian[[h[l]]](par = par, y = y)
              if(!is.null(weights))
                hessian_l <- hessian_l * weights
              adj <- adj + hessian_l * (eta[[k]] - eta_old[[k]])
            }
          }
          ew <- ew_CG[[j]]
          wj <- if(is.null(weights)) ew$weights else ew$weights * weights
          wj[!is.finite(wj)] <- 0
          wj[wj < 0] <- 0
          ew$eta <- ew$eta - adj / wj
        } else {
          ew <- .update(par = par, y = y,
            eta = if(iter[1L] > 0L) eta[[j]] else etastart[[j]],
            family = family, which = j)
        }

        ## Start inner loop.
        maxit_parameter <- if(use_CG) 1L else maxit[2L]
        while((eps[2L] > stop.eps[2L]) && (iter[2L] < maxit_parameter)) {
          ## Current log-likelihood.
          if(is.null(weights)) {
            ll0 <- log_likelihood(par = map2par(eta), y = y)
          } else {
            ll0 <- sum(pdf(par = map2par(eta), y = y, log = TRUE) * weights, na.rm = TRUE)
          }
          ll02 <- ll0
          eta_parameter <- eta
          fit_parameter <- fit[[j]]
          sfit_parameter <- sfit[[j]]
          penalty_parameter <- penalty[j]

          ## Fit linear part.
          if(length(xterms[[j]])) {
            ## Compute partial residuals.
            eta[[j]] <- eta[[j]] - fit[[j]]$fitted.values
            e <- ew$eta - eta[[j]]

            ## Weights.
            wj <- if(is.null(weights)) ew$weights else ew$weights * weights
            wj[!is.finite(wj)] <- 0
            wj[wj < 0] <- 0

            ## Design matrix.
            Xj <- x[, xterms[[j]], drop = FALSE]

            ## Estimate weighted linear model.
            if(ridge) {
              m <- ridge.lm.wfit(Xj, e, wj, penalty = penalty[j], control)
              penalty[j] <- m$penalty
            } else {
              m <- rs_lm_wfit(Xj, e, wj, linear.cache[[j]])
              if(!stepwise_candidate)
                m$vcov <- vcov_lm_wfit_safe(m)
            }

            ## If linear model does not improve the fit, use ML.
            etai <- eta
            etai[[j]] <- etai[[j]] + m$fitted.values

            ll1 <- candidate_log_likelihood(etai)

            if(ll1 < ll02 && isTRUE(control$backup)) {
              ll <- function(par) {
                etai <- eta
                etai[[j]] <- etai[[j]] + drop(Xj %*% par)
                -candidate_log_likelihood(etai) + lambda * sum(par^2)
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
                etai <- eta
                etai[[j]] <- etai[[j]] + m$fitted.values
                ll1 <- candidate_log_likelihood(etai)
              }
            }

            ## Update predictor.
            fit[[j]]$fitted.values <- m$fitted.values
            fit[[j]]$coefficients <- m$coefficients

            if(ridge)
              fit[[j]]$penalty <- m$penalty

            if(!is.null(m$edf))
              fit[[j]]$edf <- m$edf

            if(!stepwise_candidate) {
              if(!is.null(m$vcov)) {
                fit[[j]]$vcov <- m$vcov
              } else {
                Xw <- Xj * sqrt(pmax(wj, 0))
                XWX <- crossprod(Xw)

                ridge.eps <- 1e-8 * mean(diag(XWX))
                if(!is.finite(ridge.eps) || ridge.eps <= 0) ridge.eps <- 1e-8
                diag(XWX) <- diag(XWX) + ridge.eps

                R <- tryCatch(chol(XWX), error = function(e) NULL)

                if(is.null(R)) {
                  ridge.lambda <- ridge.eps
                  for(k in 1:6) {
                    Xt <- XWX + diag(ridge.lambda, ncol(XWX))
                    R <- tryCatch(chol(Xt), error = function(e) NULL)
                    if(!is.null(R)) {
                      XWX <- Xt
                      break
                    }
                    ridge.lambda <- ridge.lambda * 10
                  }
                }
                if(is.null(R)) {
                  fit[[j]]$vcov <- MASS::ginv(XWX)
                } else {
                  fit[[j]]$vcov <- chol2inv(R)
                }
              }

              colnames(fit[[j]]$vcov) <- rownames(fit[[j]]$vcov) <- colnames(Xj)
            }

            eta[[j]] <- eta[[j]] + fit[[j]]$fitted.values
          }

          ## Fit specials part.
          if(length(sterms[[j]])) {
            for(k in sterms[[j]]) {
              ## Compute partial residuals.
              eta[[j]] <- eta[[j]] - sfit[[j]][[k]]$fitted.values
              e <- ew$eta - eta[[j]]

              ## The default mgcv fitter consumes this private cache. Do not
              ## attach it to the stored special or expose it to user fitters.
              sk <- specials[[k]]
              default.fitter <- inherits(sk, "mgcv.smooth") &&
                !inherits(sk, c("smooth", "special")) &&
                is.null(sk$special.wfit)
              if(use.cache && default.fitter)
                sk$.rs_cache <- smooth.cache[[j]][[k]]

              ## Additive model term fit.
              fs <- if(is.null(weights)) {
                special.wfit(sk, e, ew$weights, y, eta, j, family, control,
                  transfer = sfit[[j]][[k]]$transfer, iter = iter)
              } else {
                special.wfit(sk, e, ew$weights * weights, y, eta, j, family, control,
                  transfer = sfit[[j]][[k]]$transfer, iter = iter)
              }

              sfit[[j]][[k]] <- fs
              sfit[[j]][[k]]$selected <- TRUE
              eta[[j]] <- eta[[j]] + fs$fitted.values
            }
          }

          ## Safeguard a complete parameter update.  Smooth and linear terms
          ## can compensate for each other and must be accepted together.
          ll1 <- candidate_log_likelihood(eta)
          if(!use_CG) {
            initial.step <- if(iter[1L] > 0L || iter[2L] > 0L)
              control$step else 1
            accepted <- safeguard(eta_parameter, eta, ll0, initial.step,
              sfit_parameter, sfit[[j]])
            if(accepted$step < 1) {
              if(accepted$accepted) {
                state <- blend_parameter(fit_parameter, fit[[j]],
                  sfit_parameter, sfit[[j]], accepted$step)
                fit[[j]] <- state$fit
                sfit[[j]] <- state$sfit
              } else {
                fit[[j]] <- fit_parameter
                sfit[[j]] <- sfit_parameter
                penalty[j] <- penalty_parameter
              }
            }
            eta <- accepted$eta
            ll1 <- accepted$logLik
            objective.accepted <- objective.accepted ||
              isTRUE(accepted$objective)
            safeguard.failed <- safeguard.failed || !accepted$accepted
            set_step(j, accepted$step)
            if(iter[1L] < 1L)
              etastart[[j]] <- eta[[j]]
          }

          ## Stopping criterion.
          eps[2L] <- abs(ll1 - ll0) / (abs(ll0) + 1e-08)

          ## Update working response.
          if((eps[2L] > stop.eps[2L]) && !use_CG) {
            par <- map2par(eta)
            ew <- .update(par = par, y = y,
              eta = if(iter[1L] > 0L) eta[[j]] else etastart[[j]],
              family = family, which = j)
          }

          ## Update inner loop iterator.
          iter[2L] <- iter[2L] + 1L
        }

        ## Reset inner iterator and stopping criterion.
        iter[2L] <- 0L
        eps[2L] <- stop.eps[2L] + 1
      }

      ## The CG correction is coupled across parameters.  Safeguard the whole
      ## sweep rather than its individual parameter updates.
      if(use_CG) {
        ll1 <- candidate_log_likelihood(eta)
        initial.step <- if(iter[1L] > 0L || iter_outer > 0L)
          control$step else 1
        accepted <- safeguard(eta_sweep, eta, outer_ll0, initial.step,
          sfit_sweep, sfit)
        if(accepted$step < 1) {
          if(accepted$accepted) {
            for(j in np) {
              if(control$fixed[[j]])
                next
              state <- blend_parameter(fit_sweep[[j]], fit[[j]],
                sfit_sweep[[j]], sfit[[j]], accepted$step)
              fit[[j]] <- state$fit
              sfit[[j]] <- state$sfit
            }
          } else {
            fit <- fit_sweep
            sfit <- sfit_sweep
            penalty <- penalty_sweep
          }
        }
        eta <- accepted$eta
        ll1 <- accepted$logLik
        objective.accepted <- objective.accepted ||
          isTRUE(accepted$objective)
        safeguard.failed <- safeguard.failed || !accepted$accepted
        set_all_steps(accepted$step)
      }

      ## For Cole and Green.
      iter_outer <- iter_outer + 1L
      if(use_CG)
        eps_outer <- abs(ll1 - outer_ll0) / (abs(outer_ll0) + 1e-08)
    }

    ## New log-likelihood.
    if(is.null(weights)) {
      llo1 <- log_likelihood(par = map2par(eta), y = y)
    } else {
      llo1 <- sum(pdf(par = map2par(eta), y = y, log = TRUE) * weights, na.rm = TRUE)
    }

    ## Stopping criterion.
    eps[1L] <- abs(llo1 - llo0) / (abs(llo0) + 1e-08)

    ## A material decrease should have been prevented by the safeguard.
    if(!objective.accepted && iter[1L] > 0L && is.finite(llo0) &&
        is.finite(llo1) && llo1 < llo0) {
      rel.inc <- (llo0 - llo1) / (abs(llo0) + 1e-08)
      if(rel.inc > stop.eps[1L]) {
        dev.warn <- dev.warn + 1L
        dev.warn.max <- max(dev.warn.max, rel.inc)
      }
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
      itxt <- paste0(paste0("GAMLSS-", if(use_CG) "CG" else "RS", " iteration "),
        fmt(iter[1L], nchar(as.character(maxit[1L])), digits = 0),
        ": Global Deviance = ", round(-2 * llo1, digits = 4),
        " eps = ", fmt(eps[1L], width = 8, digits = 8), "    ")
      cat(itxt, if(control$flush) NULL else "\n")
    }
  }

  if(control$trace & control$flush)
    cat("\n")

  if(dev.warn > 0L) {
    warning(sprintf(
      paste0("Global deviance increased materially in ",
        "%d outer iteration(s) (maximum relative increase %.3g); ",
        "check convergence or reduce the step length."),
      dev.warn, dev.warn.max
    ))
  }

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
          sfit[[j]][[i]]$.from_start <- NULL
          sfit[[j]][[i]]$.rs_smoothing <- NULL
          if(isTRUE(control$light) || stepwise_candidate) {
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
  if(length(sfit)) {
    for(j in names(sfit)) {
      for(i in names(sfit[[j]])) {
        if(!is.null(sfit[[j]][[i]]$transfer$names)) {
          if(length(sfit[[j]][[i]]$transfer$names) == length(sfit[[j]][[i]]$coefficients)) {
            names(sfit[[j]][[i]]$coefficients) <- sfit[[j]][[i]]$transfer$names
          }
        }
      }
    }
  }

  if(ridge) {
    attr(fit, "edf") <- unlist(sapply(fit, function(x) x$edf))
  }

  if(stepwise_candidate && length(fit)) {
    for(j in names(fit))
      fit[[j]]$fitted.values <- NULL
  }

  ## Message if not converged due to NAs or Inf!
  par <- map2par(eta)
  d <- pdf(par = par, y = y, log = TRUE)
  if(!is.null(weights))
    d <- d * weights
  if(any(is.na(d))) {
    warning("NA log-density values in the last iteration of the RS algorithm!")
  }
  if(any(!is.finite(d))) {
    warning("non-finite log-density values in the last iteration of the RS algorithm!")
  }

  ## A scale smooth can expose an unbounded likelihood by approaching zero at
  ## one or a few observations. Warn about the fitted function without imposing
  ## a distribution-independent lower bound on sigma.
  if(control$sigma.tol > 0 && "sigma" %in% names(par) &&
      length(sterms[["sigma"]])) {
    sigma <- par[["sigma"]]
    sigma <- sigma[is.finite(sigma) & sigma > 0]
    if(length(sigma)) {
      reference <- median(sigma)
      ratio <- min(sigma) / reference
      if(is.finite(ratio) && ratio < control$sigma.tol) {
        warning(sprintf(paste0(
          "fitted sigma approaches zero relative to its median ",
          "(minimum/median = %.3g); the likelihood may be near-singular ",
          "and smooth estimates unstable. Consider a stronger smoothing ",
          "criterion or a smaller basis dimension."), ratio))
      }
    }
  }

  converged <- is.finite(eps[1L]) && eps[1L] <= stop.eps[1L] &&
    !safeguard.failed

  rval <- list(
    "fitted.values" = if(stepwise_candidate) NULL else as.data.frame(eta),
    "fitted.specials" = sfit,
    "fitted.linear" = fit,
    "coefficients" = coef_lin,
    "iterations" = iter[1L],
    "converged" = converged,
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

## Memoize the last evaluation of an opted-in, deterministic family callback.
## Keep IEEE values and attributes distinct, including signed zero. Never
## suppress warnings by reusing an evaluation that produced one.
rs_cached_function <- function(fun)
{
  force(fun)
  last.args <- last.value <- NULL
  valid <- FALSE
  function(...) {
    args <- list(...)
    if(valid && identical(args, last.args, num.eq = FALSE, single.NA = FALSE))
      return(last.value)
    valid <<- FALSE
    warned <- FALSE
    value <- withCallingHandlers(fun(...), warning = function(w) warned <<- TRUE)
    if(!warned) {
      last.args <<- args
      last.value <<- value
      valid <<- TRUE
    }
    value
  }
}

## Reuse the same LINPACK QR and residual calculation as lm.wfit(). Restrict
## reuse to full-rank fits with strictly positive weights; all other cases
## keep lm.wfit()'s pivoting and zero-weight handling.
rs_lm_wfit <- function(x, y, w, cache = NULL)
{
  if(is.environment(cache) && !is.null(cache$fit) &&
      identical(w, cache$w, num.eq = FALSE, single.NA = FALSE) &&
      identical(x, cache$x, num.eq = FALSE, single.NA = FALSE)) {
    m <- cache$fit
    wy <- y * cache$sqrtw
    m$coefficients <- qr.coef(m$qr, wy)
    m$residuals <- qr.resid(m$qr, wy) / cache$sqrtw
    m$fitted.values <- y - m$residuals
    effects <- qr.qty(m$qr, wy)
    names(effects) <- names(m$effects)
    m$effects <- effects
    return(m)
  }

  m <- lm.wfit(x, y, w, method = "qr")
  if(is.environment(cache)) {
    cache$fit <- NULL
    if(is.null(dim(y)) && m$rank == ncol(x) && ncol(x) > 0L &&
        all(is.finite(w)) && all(w > 0)) {
      cache$x <- x
      cache$w <- w
      cache$sqrtw <- sqrt(w)
      cache$fit <- m
    }
  }
  m
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

## Function to check values of score and hessian vectors
deriv_checks <- function(x, is.weight = FALSE)
{
  ## Scores and weights are normally already finite and within bounds. Avoid
  ## several full-vector subassignments in that common case.
  if(!length(x))
    return(x)
  if(!anyNA(x)) {
    if(is.weight) {
      if(min(x) >= 1e-10 && max(x) <= 1e+10)
        return(x)
    } else if(min(x) >= -1e+10 && max(x) <= 1e+10) {
      return(x)
    }
  }

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

## Compute working response z and weights hessian from family.
.update <- function(par, y, eta, family, which)
{
  if(is.null(family$update)) {
    score <- deriv_checks(family$score[[which]](par = par, y = y, id = which), is.weight = FALSE)
    hessian <- deriv_checks(family$hessian[[which]](par = par, y = y, id = which), is.weight = TRUE)
    z <- eta + 1 / hessian * score
    return(list("eta" = z, "weights" = hessian))
  } else {
    return(family$update(par = par, y = y, eta = eta, which = which))
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

  if(!is.finite(penalty) || penalty < 0)
    penalty <- 1

  w[!is.finite(w)] <- 0
  w[w < 0] <- 0

  nc <- ncol(x)
  n <- length(w)

  i <- which(colnames(x) == "(Intercept)")
  I <- rep(1, nc)
  if(length(i))
    I[i] <- 0.0

  only_itcpt <- all(colnames(x) == "(Intercept)")

  sw <- sqrt(w)
  XW <- x * sw
  XWX <- crossprod(XW)
  XWy <- crossprod(XW, y * sw)

  fp <- function(pen, rf = FALSE) {
    pen <- max(pen, 0)

    S <- diag(I * pen, nc)
    Q <- XWX + S
    Q <- 0.5 * (Q + t(Q))

    eps <- 1e-8 * mean(diag(Q))
    if(!is.finite(eps) || eps <= 0)
      eps <- 1e-8
    diag(Q) <- diag(Q) + eps

    cholQ <- tryCatch(chol(Q), error = function(e) NULL)

    if(is.null(cholQ)) {
      lam <- eps
      for(k in 1:10) {
        Qt <- Q + diag(lam, nc)
        cholQ <- tryCatch(chol(Qt), error = function(e) NULL)
        if(!is.null(cholQ)) {
          Q <- Qt
          break
        }
        lam <- lam * 10
      }
    }

    if(is.null(cholQ)) {
      ev <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
      shift <- max(1e-8, -min(ev) + 1e-6)
      Q <- Q + diag(shift, nc)
      cholQ <- chol(Q)
    }

    b <- backsolve(cholQ, forwardsolve(t(cholQ), XWy))
    b <- drop(b)

    fit <- drop(x %*% b)

    Tmat <- backsolve(cholQ, forwardsolve(t(cholQ), XWX))
    edf <- sum(diag(Tmat))

    if(!is.finite(edf))
      edf <- nc
    if(edf > (n - 1e-8))
      edf <- n - 1e-8

    if(rf) {
      names(b) <- colnames(x)
      return(list(
        coefficients = b,
        fitted.values = fit,
        edf = edf,
        penalty = pen,
        vcov = if(isTRUE(control$.stepwise_candidate)) NULL else chol2inv(cholQ),
        df = n - edf
      ))
    }

    rss <- sum(w * (y - fit)^2)

    switch(tolower(control$criterion),
      "gcv"  = rss * n / (n - edf)^2,
      "aic"  = rss + 2 * edf,
      "gaic" = rss + K * edf,
      "aicc" = rss + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1),
      "bic"  = rss + log(n) * edf
    )
  }

  if(!only_itcpt) {
    lower <- max(1e-8, penalty / 10)
    upper <- max(lower * 10, penalty * 10, 1e-6)
    opt <- nlminb(penalty, objective = fp, lower = lower, upper = upper)
  } else {
    opt <- list(par = 0.0)
  }

  fp(opt$par, rf = TRUE)
}

## Safe vcov.
vcov_lm_wfit_safe <- function(m, ridge = 1e-8, maxit = 6,
  ginv_fallback = TRUE, warn_negative_weights = TRUE)
{
  if(!is.null(m$weights)) {
    if(any(m$weights < 0, na.rm = TRUE)) {
      msg <- "Negative weights encountered."
      if(warn_negative_weights) warning(msg) else stop(msg)
    }
  }

  r <- m$rank
  p <- length(m$coefficients)
  piv <- m$qr$pivot

  if(r == 0L) {
    V <- matrix(NA_real_, p, p)
    return(V)
  }

  R <- qr.R(m$qr)[1:r, 1:r, drop = FALSE]

  ## Form crossprod(R) = X'WX on the estimable subspace.
  XtWX <- crossprod(R)

  eps <- ridge * mean(diag(XtWX))
  if(!is.finite(eps) || eps <= 0) eps <- ridge
  diag(XtWX) <- diag(XtWX) + eps

  Rc <- tryCatch(chol(XtWX), error = function(e) NULL)

  if(is.null(Rc)) {
    lambda <- eps
    for(k in seq_len(maxit)) {
      Xt <- XtWX + diag(lambda, ncol(XtWX))
      Rc <- tryCatch(chol(Xt), error = function(e) NULL)
      if(!is.null(Rc)) {
        XtWX <- Xt
        break
      }
      lambda <- lambda * 10
    }
  }

  Vinv <- if(is.null(Rc)) {
    if(ginv_fallback) MASS::ginv(XtWX) else stop("Cannot invert weighted crossproduct.")
  } else {
    chol2inv(Rc)
  }

  V <- matrix(NA_real_, p, p)
  V[piv[1:r], piv[1:r]] <- Vinv
  V
}
