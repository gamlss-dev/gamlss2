## Function to compute multivariate normal samples
## given the maximum likelihood estimator.
sampling <- function(object, R = 100, ...)
{
  V <- vcov(object, ...)
  if(any(eigen(V)$values < 0)) {
    V <- as.matrix(Matrix::nearPD(V)$mat)
  }
  Cv <- chol(V)
  cb <- coef(object, ...)
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
  iterthin <- as.integer(seq(burnin, n.iter, by = thinning))
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
          if(!inherits(specials[[i]], "mgcv.smooth"))
            stop("only mgcv smooth terms are allowed!")
          sfit[[j]][[i]] <- list("fitted.values" = rep(0.0, n), "edf" = 0.0, "selected" = FALSE)
          samples[[j]]$s[[i]] <- matrix(NA, nrow = nsave,
            ncol = ncol(specials[[i]]$X) + length(specials[[i]]$S) + 1L)
          if(!is.null(cstart)) {
            sj <- grep(paste0(j, ".s.", i), names(cstart), fixed = TRUE, value = TRUE)
            if(length(sj)) {
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
    if(nes[[j]])
      etastart[[j]] <- eta[[j]]
  }

  ## For printing.
  if(control$flush) {
    control$flush <- interactive()
  }

  ## Start MCMC.
  if(control$trace) {
    if(!is.null(control$light)) {
      if(control$light)
        cat("Start sampling ...\n")
    }
  }
  
  for(iter in 1:n.iter) {
    for(j in np) {
      ## Check if paramater is fixed.
      if(control$fixed[[j]])
        stop("fixed parameters not supported yet!")

      ## Outer loop working response and weights.
      peta <- if(iter > 1L) {
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
    }
  }
}

