## Function simply evaluates the
## special terms in the model formula
## and assigns appropriate model fitting
## functions for the backfitting steps.
special_terms <- function(x, data, binning = FALSE, digits = Inf, ...)
{
  sterms <- list()

  if(length(x)) {
    for(j in unlist(x)) {
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
        smooth.construct.wfit
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

## Fitting function for mgcv smooth terms.
smooth.construct.wfit <- function(x, z, w, y, eta, j, family, control, transfer, iter)
{
  ## Number of observations.
  n <- length(z)

  if(control$binning) {
    rw <- numeric(length(x$binning$nodups))
    rz <- numeric(length(x$binning$nodups))
    calc_Xe(x$binning$sorted.index, w, z, rw, rz, x$binning$order)
  }

  ## Pre compute matrices.
  if(control$binning) {
    XWz <- crossprod(x$X, rz)
    XWX <- calc_XWX(x$X, 1/rw, x$sparse_index)
  } else {
    XW <- x$X * w
    XWX <- crossprod(XW, x$X)
    XWz <- crossprod(XW, z)
  }
  S <- diag(1e-05, ncol(x$X))

  if(!is.null(x$control)) {
    control[names(x$control)] <- x$control
    if(!is.null(control$method))
      control$criterion <- tolower(control$method)
  }
  if(is.null(control$criterion))
    control$criterion <- "aicc"

  ## Set up smoothing parameters.
  if(iter[1L] > -1) {
    lambdas <- transfer$lambdas
  } else {
    lambdas <- 10
  }
  if(is.null(lambdas)) {
    lambdas <- if(is.null(control$start)) 10 else control$start
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

  if(control$criterion == "ml" & length(x$S) < 2L & !control$binning & localML) {
    ## Local ML method, only for pb2() yet!
    order <- x$m[1L]
    if(is.null(order))
      order <- 1

    N <- sum(w != 0)

    for(it in 1:50) {
      P <- try(chol2inv(chol(XWX + lambdas * x$S[[1L]])), silent = TRUE)
      if(inherits(P, "try-error"))
        P <- solve(XWX + lambdas * x$S[[1L]])

      b <- drop(P %*% XWz)
      fit <- drop(x$X %*% b)

      if(control$binning)
        fit <- fit[x$binning$match.index]

      edf <- sum(diag(XWX %*% P))

      sig2 <- sum(w * (z - fit)^2) / (N - edf)
      tau2 <- drop(t(b) %*% x$S[[1L]] %*% b) / (edf - order)

      if(tau2 < 1e-07) tau2 <- 1e-07
      lambdas.old <- lambdas
      lambdas <- sig2/tau2
      if(lambdas < 1e-07) lambdas <- 1e-07
      if(lambdas > 1e+07) lambdas <- 1e+07
      if(abs(lambdas - lambdas.old) < 1e-07 || lambdas > 1e+10) break
    }

    return(list("coefficients" = b, "fitted.values" = fit, "edf" = edf,
      "lambdas" = lambdas, "vcov" = P, "df" = n - edf))
  } else {
    ## Function to search for smoothing parameters using GCV etc.
    fl <- function(l, rf = FALSE) {
      for(j in 1:length(x$S))
        S <- S + l[j] * x$S[[j]]

      P <- try(chol2inv(chol(XWX + S)), silent = TRUE)
      if(inherits(P, "try-error"))
        P <- solve(XWX + S)

      b <- drop(P %*% XWz)

      fit <- drop(x$X %*% b)

      if(control$binning)
        fit <- fit[x$binning$match.index]

      edf <- sum(diag(XWX %*% P))

      if(rf) {
        return(list("coefficients" = b, "fitted.values" = fit, "edf" = edf,
          "lambdas" = l, "vcov" = P, "df" = n - edf))
      } else {
        if(isTRUE(control$logLik)) {
          eta[[j]] <- eta[[j]] + fit
          rss <- family$loglik(y, family$map2par(eta))
        } else {
          rss <- sum(w * (z - fit)^2)
        }

        rval <- switch(tolower(control$criterion),
          "gcv" = rss * n / (n - edf)^2,
          "aic" = rss + 2 * edf,
          "gaic" = rss + K * edf,
          "aicc" = rss + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1),
          "bic" = rss + log(n) * edf
        )

        return(rval)
      }
    }

    if(is.null(x$sp)) {
      opt <- nlminb(lambdas, objective = fl, lower = lambdas / 10, upper = lambdas * 10)
    } else {
      opt <- list(par = x$sp)
    }
  }

  rval <- fl(opt$par, rf = TRUE)

  rval$transfer <- list("lambdas" = rval$lambdas)

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
    if(is.null(x$model))
      return(predict(x, newdata = data))
    else
      return(predict(x$model, newdata = data))
  }
}

