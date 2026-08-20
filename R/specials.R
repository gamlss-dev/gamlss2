## Function simply evaluates the
## special terms in the model formula
## and assigns appropriate model fitting
## functions for the backfitting steps.
special_terms <- function(x, data, binning = FALSE, digits = Inf, ...)
{
  sterms <- list()
  if(length(x)) {
    for(j in unique(unlist(x))) {
      vj <- all.vars(parse(text = j))
      vj <- vj[vj %in% names(data)]

      if(length(vj) < 1L) {
        vjf <- as.formula(paste0("~", j))
        vjf <- fake_formula(vjf)
        vj <- attr(terms(vjf), "term.labels")
      }

      binj <- binning

      if(binj) {
        dj <- data[, vj, drop = FALSE]
        if(is.finite(digits)) {
          for(v in vj) {
            if(is.numeric(dj[[v]])) {
              dj[[v]] <- round(dj[[v]], digits = digits)
            }
          }
        }
        dj <- apply(dj, 1, paste, sep = "\r", collapse = ";")

        bn <- list()
        bn$nodups <- which(!duplicated(dj))
        bn$match.index <- match(dj, dj[bn$nodups])
        bn$order <- order(bn$match.index)
        bn$sorted.index <- bn$match.index[bn$order]

        dj <- data[bn$nodups, vj, drop = FALSE]
      }

      ## Change constructor if possible.
      sjp <- eval(parse(text = paste0("quote(", j, ")")))
      sjpc <- as.character(sjp[1L])
      changed <- FALSE
      if(sjpc == "pb" & FALSE) {
        sjp[[1L]] <- as.name("pb2")
        j <- deparse(sjp)
        changed <- TRUE
      }

      sj <- eval(parse(text = j), envir = if(binj) dj else data)

      ## For class "smooth", binning is not possible.
      if(inherits(sj, "smooth") & binning) {
        warning(paste0("binning is not possible for 'smooth' term ", j, "!"))
        sj <- eval(parse(text = j), envir = data)
        binj <- FALSE
      }

      if(any(grepl(".smooth.spec", class(sj)))) {
        stopifnot(requireNamespace("mgcv"))
        knots <- list(...)$knots

        absorb.cons <- if(is.null(sj$xt$absorb.cons)) TRUE else isTRUE(sj$xt$absorb.cons)
        scale.penalty <- if(is.null(sj$xt$scale.penalty)) TRUE else isTRUE(sj$xt$scale.penalty)

        select <- isTRUE(list(...)$select)

        sj <- mgcv::smoothCon(sj, data = if(binj) dj else data, knots = knots,
          absorb.cons = absorb.cons, scale.penalty = scale.penalty,
          null.space.penalty = select)

        for(i in 1:length(sj)) {
          sj[[i]]$orig.label <- j
          if(binj) {
            sj[[i]]$binning <- bn
            sj[[i]]$sparse_index <- calc_sparse_index(sj[[i]]$X)
          }
          if(select) {
            sj[[i]]$S[[length(sj[[i]]$S)]] <- sj[[i]]$S[[length(sj[[i]]$S)]] + diag(1/sqrt(ncol(sj[[i]]$X)), ncol(sj[[i]]$X))
          }
        }
        sjn <- sapply(sj, function(x) x$label)
        if(changed) {
          sjn <- j
          sj[[1L]]$label <- sjn
          sj[[1L]]$orig.label <- sjn
        }
        names(sj) <- sjn
        sterms <- c(sterms, sj)
      } else {
        if(binj) {
          sj$binning <- bn
        }
        if(inherits(sj, "matrix")) {
          sj <- list("X" = sj, "S" = list(diag(1, ncol(sj))),
            "label" = j, "term" = all.vars(parse(text = j)), "dim" = 1L)
        }
        sterms[[j]] <- sj
      }
    }
  }

  if(any(dups <- duplicated(names(sterms)))) {
    dups <- which(dups)
    for(j in seq_along(dups)) {
      dn <- names(sterms)[dups[j]]
      dn <- paste0(dn, ".", j)
      names(sterms)[dups[j]] <- dn
    }
  }

  return(sterms)
}

## Calculate index matrix of non-zero elements.
calc_sparse_index <- function(x, ...)
{
  if(is.null(dim(x)))
    return(NULL)
  index <- apply(x, 1, function(x) {
    which(x != 0)
  })
  if(length(index) < 1)
    return(NULL)
  if(is.list(index)) {
    n <- max(sapply(index, length))
    index <- lapply(index, function(x) {
      if((nx <- length(x)) < n)
        x <- c(x, rep(-1L, length = n - nx))
      x
    })
    index <- do.call("rbind", index)
  } else {
    index <- if(is.null(dim(index))) {
      matrix(index, ncol = 1)
    } else t(index)
  }
  storage.mode(index) <- "integer"
  index
}

## Special term fit function, works with gamlss model terms, too.
special.wfit <- function(x, z, w, y, eta, j, family, control, ...)
{
  if(inherits(x, "smooth")) {
    call <- attr(x, "call")
    call[[2]] <- quote(x)
    fe <- eval(call)
    if(!is.null(fe$y))
      fe$fitted.values <- fe$y
    fit <- list(
      "fitted.values" = as.numeric(fe$fitted.values),
      "coefficients" = fe$coefSmo,
      "lambdas" = fe$lambda,
      "edf" = fe$nl.df,
      "df" = length(z) - fe$nl.df,
      "model" = fe$model
    )
  } else {
    if(inherits(x, "special")) {
      fit <- special_fit(x = x, z = z, w = w, y = y, eta = eta, j = j, family = family, control = control, ...)
    } else {
      ff <- if(is.null(x$special.wfit)) {
        smooth.construct_wfit
      } else {
        x$special.wfit
      }
      fit <- ff(x, z, w, y, eta, j, family, control, ...)
    }
  }

  return(fit)
}

## Reduced working response and weights binning.
calc_Xe <- function(ind, weights, response, rweights, rresponse, oind, uind = NULL)
{
  .Call("calc_Xe", as.integer(ind), as.numeric(weights), 
    as.numeric(response), as.numeric(rweights), as.numeric(rresponse),
    as.integer(oind), PACKAGE = "gamlss2")
}

## Fast block diagonal crossproduct with weights.
calc_XWX <- function(x, w, index = NULL)
{
  if(is.null(index)) {
    rval <- crossprod(x / w, x)
  } else {
    if(is.null(dim(index)))
      index <- matrix(index, ncol = 1)
    rval <- .Call("calc_XWX", x, w, index, PACKAGE = "gamlss2")
  }
  rval
}

## Fused dense weighted crossproducts.
calc_XWXz <- function(x, w, z)
{
  .Call("calc_XWXz", x, as.numeric(w), as.numeric(z), PACKAGE = "gamlss2")
}

## Weighted Demmler-Reinsch reparameterization for a single penalty.
smooth.construct_dr <- function(XWX, XWz, S, ridge = 1e-05)
{
  p <- ncol(XWX)
  if(!is.matrix(S) || !identical(dim(S), c(p, p)))
    stop("invalid penalty matrix")
  if(!isSymmetric(S, tol = sqrt(.Machine$double.eps)))
    stop("penalty matrix is not symmetric")

  ## The small ridge penalty is part of the metric so that the
  ## reparameterized fit is equivalent to the direct penalized solve.
  G <- XWX
  diag(G) <- diag(G) + ridge
  R <- chol(G)
  Ri <- backsolve(R, diag(p))

  S <- (S + t(S)) / 2
  St <- crossprod(Ri, S %*% Ri)
  St <- (St + t(St)) / 2
  ev <- eigen(St, symmetric = TRUE)

  tol <- sqrt(.Machine$double.eps) * max(1, abs(ev$values))
  if(any(ev$values < -tol))
    stop("penalty matrix is not positive semi-definite")
  ev$values[ev$values < 0] <- 0

  T <- Ri %*% ev$vectors
  h <- 1 - ridge * colSums(T^2)

  rval <- list(
    "T" = T,
    "d" = ev$values,
    "h" = h,
    "g" = drop(crossprod(T, XWz))
  )
  if(!all(vapply(rval, function(x) all(is.finite(x)), logical(1L))))
    stop("non-finite Demmler-Reinsch reparameterization")
  rval
}

## Fitting function for mgcv smooth terms.
smooth.construct_wfit <- function(x, z, w, y, eta, j, family, control, transfer, iter)
{
  ## Number of observations.
  n <- length(z)

  if(control$binning) {
    rw <- numeric(length(x$binning$nodups))
    rz <- numeric(length(x$binning$nodups))
    calc_Xe(x$binning$sorted.index, w, z, rw, rz, x$binning$order)
  }

  ## Pre compute matrices.
  zWz <- NULL
  if(control$binning) {
    XWz <- crossprod(x$X, rz)
    XWX <- calc_XWX(x$X, 1/rw, x$sparse_index)
  } else {
    use.symmetric <- ncol(x$X) >= 40L &&
      all(is.finite(w)) && all(w >= 0)
    if(use.symmetric) {
      crossproducts <- calc_XWXz(x$X, w, z)
      XWX <- crossproducts$XWX
      XWz <- crossproducts$XWz
      zWz <- crossproducts$zWz
    } else {
      ## Preserve the previous behavior for non-finite or negative weights.
      XW <- x$X * w
      XWX <- crossprod(XW, x$X)
      XWz <- crossprod(XW, z)
    }
  }
  S <- diag(1e-05, ncol(x$X))

  if(!is.null(x$control)) {
    control[names(x$control)] <- x$control
    if(!is.null(control$method))
      control$criterion <- tolower(control$method)
  }
  if(is.null(control$criterion))
    control$criterion <- "aicc"

  ## Extra penalty for selection.
  if(isTRUE(control$termselect)) {
    df <- ncol(x$X)
    bml <- drop(solve(XWX + diag(1e-08, df), XWz))

    pen <- function(b) {
      A <- 1 / rep(sqrt(sum(b^2)), df) * 1 / rep(sqrt(sum(bml^2)), df)
      A <- if(length(A) < 2L) matrix(A, 1, 1) else diag(A)
      A
    }

    b0 <- if(is.null(transfer$coefficients)) bml else  transfer$coefficients

    x$S[[length(x$S) + 1L]] <- pen(b0)
  }

  ## Set up smoothing parameters.
  if(iter[1L] > -1) {
    lambdas <- transfer$lambdas
  } else {
    lambdas <- 1.0
  }
  if(is.null(lambdas)) {
    lambdas <- if(is.null(control$start)) 1.0 else control$start
  }
  lambdas <- rep(lambdas, length.out = length(x$S))

  ## Penalty for AIC.
  K <- if(is.null(control$K)) 2 else control$K

  ## Local ML check.
  localML <- isTRUE(x$localML)
  if(!localML) {
    if(control$criterion == "ml")
      control$criterion <- "aicc"
  }

  ## Choose the direct or Demmler-Reinsch solver. In automatic mode the
  ## reparameterization is delayed until the current search needs a restart.
  dr.mode <- control$demmler.reinsch
  if(is.null(dr.mode))
    dr.mode <- "auto"
  if(is.character(dr.mode) && length(dr.mode) == 1L)
    dr.mode <- tolower(dr.mode)
  if(!(identical(dr.mode, "auto") || identical(dr.mode, TRUE) ||
      identical(dr.mode, FALSE)))
    stop("'demmler.reinsch' must be one of 'auto', TRUE or FALSE")

  dr.eligible <- isTRUE(x$dim == 1L) &&
      length(x$S) == 1L &&
      is.null(x$sp) &&
      !isTRUE(x$fixed) &&
      (iter[1L] > 0L || localML)

  dr_setup <- function() {
    dr.ridge <- if(control$criterion == "ml" && localML) 0 else 1e-05
    rval <- try(smooth.construct_dr(
      XWX = XWX, XWz = XWz, S = x$S[[1L]], ridge = dr.ridge
    ), silent = TRUE)
    if(inherits(rval, "try-error"))
      NULL
    else
      rval
  }

  dr <- NULL
  dr.attempted <- FALSE
  if(dr.eligible && (identical(dr.mode, TRUE) ||
      (localML && !identical(dr.mode, FALSE)))) {
    dr <- dr_setup()
    dr.attempted <- TRUE
  }

  dr_fit <- function(lambda) {
    denominator <- 1 + as.numeric(lambda[1L]) * dr$d
    if(any(!is.finite(denominator)) || any(denominator <= 0))
      stop("invalid smoothing parameter in Demmler-Reinsch fit")

    alpha <- dr$g / denominator
    list(
      "coefficients" = drop(dr$T %*% alpha),
      "edf" = sum(dr$h / denominator)
    )
  }

  if(control$criterion == "ml" & (length(x$S) < 2L) & localML) {
    ## Local ML method, only for pb2() yet!
    order <- x$m[1L]
    if(is.null(order))
      order <- 1

    N <- sum(w != 0)

    for(it in 1:50) {
      if(is.null(dr)) {
        P <- try(chol2inv(chol(XWX + lambdas * x$S[[1L]])), silent = TRUE)
        if(inherits(P, "try-error"))
          P <- solve(XWX + lambdas * x$S[[1L]])

        b <- drop(P %*% XWz)
        edf <- sum(diag(XWX %*% P))
      } else {
        drs <- dr_fit(lambdas)
        b <- drs$coefficients
        edf <- drs$edf
      }

      fit <- drop(x$X %*% b)

      if(control$binning)
        fit <- fit[x$binning$match.index]

      sig2 <- sum(w * (z - fit)^2) / (N - edf)
      tau2 <- drop(t(b) %*% x$S[[1L]] %*% b) / (edf - order)

      if(tau2 < 1e-07) tau2 <- 1e-07
      lambdas.old <- lambdas
      lambdas <- sig2/tau2
      if(lambdas < 1e-07) lambdas <- 1e-07
      if(lambdas > 1e+07) lambdas <- 1e+07
      if(abs(lambdas - lambdas.old) < 1e-07 || lambdas > 1e+10) break
    }

    if(!is.null(dr)) {
      P <- try(chol2inv(chol(XWX + lambdas.old * x$S[[1L]])), silent = TRUE)
      if(inherits(P, "try-error"))
        P <- solve(XWX + lambdas.old * x$S[[1L]])

      b <- drop(P %*% XWz)
      fit <- drop(x$X %*% b)
      if(control$binning)
        fit <- fit[x$binning$match.index]
      edf <- sum(diag(XWX %*% P))
    }

    return(list("coefficients" = b, "fitted.values" = fit, "edf" = edf,
      "lambdas" = lambdas, "vcov" = P, "df" = n - edf))
  } else {
    if(is.null(zWz))
      zWz <- sum(w * z^2)
    criterion.evaluations <- 0L

    ## Function to search for smoothing parameters using GCV etc.
    fl <- function(l, rf = FALSE) {
      if(!rf)
        criterion.evaluations <<- criterion.evaluations + 1L
      if(!is.null(dr) && !rf) {
        drs <- dr_fit(l)
        edf <- drs$edf
        b <- drs$coefficients
      } else {
        Sl <- S
        if(length(x$S)) {
          for(k in seq_along(x$S))
            Sl <- Sl + l[k] * x$S[[k]]
        }

        A <- XWX + Sl
        R <- chol(A)

        b <- drop(backsolve(
          R,
          forwardsolve(t(R), XWz)
        ))

        P <- chol2inv(R)
        edf <- sum(XWX * P)
      }

      ## Fitted values are only needed for the final result or
      ## if the full log-likelihood is used as fitting criterion.
      if(rf || isTRUE(control$logLik)) {
        fit <- drop(x$X %*% b)

        if(control$binning)
          fit <- fit[x$binning$match.index]
      }

      if(rf) {
        return(list(
          "coefficients" = b,
          "fitted.values" = fit,
          "edf" = edf,
          "lambdas" = l,
          "vcov" = P,
          "df" = n - edf
        ))
      } else {
        if(isTRUE(control$logLik)) {
          etai <- eta
          etai[[j]] <- etai[[j]] + fit
          rss <- -2 * family$log_likelihood(
            par = family$map2par(etai),
            y = y
          )
        } else {
          ## Weighted RSS from precomputed crossproducts:
          ## z'Wz - 2 b'X'Wz + b'X'WXb.
          rss <- zWz -
            2 * sum(b * XWz) +
            sum(b * (XWX %*% b))
        }

        rval <- switch(
          tolower(control$criterion),
          "gcv" = rss * n / (n - edf)^2,
          "aic" = rss + 2 * edf,
          "gaic" = rss + K * edf,
          "aicc" = rss + 2 * edf +
            (2 * edf * (edf + 1)) / (n - edf - 1),
          "bic" = rss + log(n) * edf
        )

        return(rval)
      }
    }

    ## Check for fx = TRUE.
    if(isTRUE(x$fixed)) {
      if(is.null(x$sp)) {
        np <- if(length(x$S)) length(x$S) else 1L
        x$sp <- rep(1e-10, np)
      }
    }

    if(is.null(x$sp)) {
      rho <- log(pmax(lambdas, 1e-10))
      eps <- Inf
      lk <- 0L
      while((eps > 0.000001) && (lk < 5L)) {
        rho0 <- rho
        opt <- nlminb(
          rho,
          objective = function(rho) fl(exp(rho)),
          lower = pmax(rho - log(10), log(1e-10)),
          upper = pmin(rho + log(10), log(1e+10))
        )
        rho <- opt$par
        eps <- max(abs(rho - rho0))
        lk <- lk + 1L

        ## Reaching the search-window boundary predicts another costly restart.
        if(identical(dr.mode, "auto") && dr.eligible && !dr.attempted &&
            lk < 5L && eps >= 0.9 * log(10)) {
          dr <- dr_setup()
          dr.attempted <- TRUE
        }
      }
      opt <- list(par = exp(rho))
    } else {
      opt <- list(par = x$sp)
    }
  }

  rval <- fl(opt$par, rf = TRUE)

  rval$transfer <- list(
    "lambdas" = rval$lambdas,
    "coefficients" = rval$coefficients,
    "criterion.evaluations" = criterion.evaluations,
    "demmler.reinsch" = !is.null(dr)
  )

  return(rval)
}

## A method for fitting special terms.
special_fit <- function(x, ...)
{
  UseMethod("special_fit")
}

## A method for predicting special terms.
special_predict <- function(x, ...)
{
  UseMethod("special_predict")
}

## Default method.
special_predict.default <- function(x, data, ...)
{
  if(is.null(x)) {
    return(rep(0, nrow(data)))
  } else {
    if(is.null(x$model)) {
      return(predict(x, newdata = data))
    } else {
      return(predict(x$model, newdata = data))
    }
  }
}

## Specials extractor function after fitting the model.
specials <- function(object, model = NULL, terms = NULL, elements = NULL, ...)
{
  if(is.null(object$fitted.specials)) {
    return(NULL)
  }

  ## Extract response name, sometimes needed.
  rn <- response_name(object)

  ## Which parameter model to predict?
  if(is.null(model)) {
    model <- list(...)$what
    if(is.null(model))
      model <- list(...)$parameter
    if(is.null(model))
      model <- object$family$names
  }
  if(!is.character(model))
    model <- object$family$names[model]
  model <- object$family$names[pmatch(model, object$family$names)]

  rval <- NULL

  for(i in model) {
    if(!is.null(object$fitted.specials[[i]])) {
      it <- if(is.null(terms)) {
        names(object$fitted.specials[[i]])
      } else {
        grep2(terms, names(object$fitted.specials[[i]]), value = TRUE, fixed = TRUE)
      }

      tmp <- object$fitted.specials[[i]][it]
      names(tmp) <- paste0(i, ".", it)

      if(!is.null(elements)) {
        for(j in seq_along(tmp)) {
          cj <- class(tmp[[j]])
          if(!is.null(elements)) {
            if((length(elements) == 1L) && (elements == "names")) {
              tmp[[j]] <- names(tmp[[j]])
            } else {
              wj <- grep2(elements, names(tmp[[j]]), ignore.case = FALSE, value = TRUE, fixed = TRUE)
              if(length(wj)) {
                tmp[[j]] <- if(length(wj) > 1L) tmp[[j]][wj] else tmp[[j]][[wj]]
              }
            }
          }
        }
      }

      for(j in seq_along(tmp)) {
        if(is.list(tmp[[j]])) {
          if(is.null(tmp[[j]]$response_name)) {
            tmp[[j]]$response_name <- rn
          }
        }
      }

      rval <- c(rval, tmp)
    }
  }

  drop <- list(...)$drop
  if(is.null(drop))
    drop <- TRUE

  if((length(rval) < 2L) && drop)
    rval <- rval[[1L]]

  return(rval)
}

