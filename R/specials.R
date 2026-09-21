## Function simply evaluates the
## special terms in the model formula
## and assigns appropriate model fitting
## functions for the backfitting steps.
special_terms <- function(x, data, binning = FALSE, digits = Inf, ...)
{
  sterms <- list()
  dots <- list(...)
  if(length(x)) {
    for(j in unique(unlist(x))) {
      expr <- str2lang(j)
      vj <- all.vars(expr)
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
        ## apply(..., 1, paste) constructs and invokes an R closure once per
        ## row. Preserve its as.matrix() coercion and row-key representation,
        ## but paste whole columns in one vectorized call.
        dj <- as.matrix(dj)
        if(ncol(dj) == 1L) {
          dj <- as.character(dj[, 1L])
        } else {
          dj <- do.call(paste, c(
            lapply(seq_len(ncol(dj)), function(i) dj[, i]),
            list(sep = ";")
          ))
        }

        bn <- list()
        bn$nodups <- which(!duplicated(dj))
        bn$match.index <- match(dj, dj[bn$nodups])
        bn$order <- order(bn$match.index)
        bn$sorted.index <- bn$match.index[bn$order]

        dj <- data[bn$nodups, vj, drop = FALSE]
      }

      ## Change constructor if possible.
      sjp <- expr
      sjpc <- as.character(sjp[1L])
      changed <- FALSE
      if(sjpc == "pb" & FALSE) {
        sjp[[1L]] <- as.name("pb2")
        j <- deparse(sjp)
        changed <- TRUE
      }

      sj <- eval(sjp, envir = if(binj) dj else data)

      ## For class "smooth", binning is not possible.
      if(inherits(sj, "smooth") & binning) {
        warning(paste0("binning is not possible for 'smooth' term ", j, "!"))
        sj <- eval(sjp, envir = data)
        binj <- FALSE
      }

      if(any(grepl(".smooth.spec", class(sj)))) {
        stopifnot(requireNamespace("mgcv"))
        knots <- dots$knots

        absorb.cons <- if(is.null(sj$xt$absorb.cons)) TRUE else isTRUE(sj$xt$absorb.cons)
        scale.penalty <- if(is.null(sj$xt$scale.penalty)) TRUE else isTRUE(sj$xt$scale.penalty)

        select <- isTRUE(dots$select)

        sj <- mgcv::smoothCon(sj, data = if(binj) dj else data, knots = knots,
          absorb.cons = absorb.cons, scale.penalty = scale.penalty,
          null.space.penalty = select)

        for(i in seq_along(sj)) {
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
            "label" = j, "term" = all.vars(expr), "dim" = 1L)
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
  .Call(C_calc_Xe, as.integer(ind), as.numeric(weights),
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
    rval <- .Call(C_calc_XWX, x, w, index, PACKAGE = "gamlss2")
  }
  rval
}

## Fused dense weighted crossproducts.
calc_XWXz <- function(x, w, z, XWX = NULL)
{
  if(!is.null(XWX))
    return(.Call(C_calc_XWXz_cached, x, as.numeric(w), as.numeric(z), XWX,
      PACKAGE = "gamlss2"))
  .Call(C_calc_XWXz, x, as.numeric(w), as.numeric(z), PACKAGE = "gamlss2")
}

## Fused direct smooth-fit criterion and final-state kernel.
calc_smooth_wfit <- function(XWX, XWz, penalties, lambda, ridge,
  zWz, n, K, criterion, final = FALSE, gradient = FALSE, state = FALSE)
{
  ## mgcv represents an unpenalized (fx = TRUE) smooth with S = NULL.
  if(is.null(penalties))
    penalties <- list()
  .Call(
    C_calc_smooth_wfit,
    XWX, as.numeric(XWz), penalties, as.numeric(lambda),
    as.numeric(ridge), as.numeric(zWz), as.numeric(n), as.numeric(K),
    as.integer(criterion),
    if(isTRUE(gradient)) { if(isTRUE(state)) 4L else 2L } else
      if(isTRUE(state)) 3L else as.logical(final),
    PACKAGE = "gamlss2"
  )
}

## Native gradient from a previously computed direct fit.
calc_smooth_wfit_gradient <- function(XWX, XWz, penalties, lambda, state,
  ridge, n, K, criterion, root = NULL)
{
  .Call(C_calc_smooth_wfit_gradient_root, XWX, as.numeric(XWz), penalties,
    as.numeric(lambda), state, as.numeric(ridge), as.numeric(n),
    as.numeric(K), as.integer(criterion), root, PACKAGE = "gamlss2")
}

## Observation-space residuals for cancellation-sensitive criteria.
calc_smooth_residual <- function(X, b, z, w, index = NULL, b1 = NULL, native = TRUE)
{
  if(isTRUE(native) && is.double(X))
    return(.Call(C_calc_smooth_residual, X, as.numeric(b), as.numeric(z),
      as.numeric(w), index, if(is.null(b1)) NULL else as.numeric(b1),
      PACKAGE = "gamlss2"))
  fit <- drop(X %*% b)
  if(!is.null(index)) {
    fit <- fit[index]
    X <- X[index, , drop = FALSE]
  }
  residual <- fit - z
  rp <- drop(crossprod(X, w * residual))
  list(rss = sum(w * residual^2), residual = rp,
    rss1 = if(is.null(b1)) NULL else 2 * sum(rp * b1))
}

## Weighted Demmler-Reinsch reparameterization for a single penalty.
smooth.construct_dr <- function(XWX, XWz, S, ridge = 1e-05,
  rank = NULL, native = TRUE)
{
  p <- ncol(XWX)
  if(!is.matrix(S) || !identical(dim(S), c(p, p)))
    stop("invalid penalty matrix")
  if(!isSymmetric(S, tol = sqrt(.Machine$double.eps)))
    stop("penalty matrix is not symmetric")

  if(!is.null(rank) && (length(rank) != 1L || !is.finite(rank) ||
      rank < 0 || rank > p || rank != as.integer(rank)))
    stop("invalid penalty rank")
  if(isTRUE(native) && is.double(XWX) && is.double(S)) {
    rval <- .Call(C_calc_smooth_dr, XWX, as.numeric(XWz), S,
      as.numeric(ridge), if(is.null(rank)) -1L else as.integer(rank),
      PACKAGE = "gamlss2")
    if(!all(vapply(rval, function(x) all(is.finite(x)), logical(1L))))
      stop("non-finite Demmler-Reinsch reparameterization")
    return(rval)
  }

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
  if(!is.null(rank) && rank < p) {
    null <- seq.int(rank + 1L, p)
    if(any(abs(ev$values[null]) > tol))
      stop("penalty rank does not match its null space")
    ev$values[null] <- 0
  }

  T <- Ri %*% ev$vectors
  h <- 1 - ridge * colSums(T^2)
  for(j in which(h < sqrt(.Machine$double.eps)))
    h[j] <- sum(T[, j] * (XWX %*% T[, j]))

  rval <- list(
    "T" = T,
    "d" = ev$values,
    "h" = h,
    "M" = crossprod(T),
    "ridge" = ridge,
    "g" = drop(crossprod(T, XWz))
  )
  if(!all(vapply(rval, function(x) all(is.finite(x)), logical(1L))))
    stop("non-finite Demmler-Reinsch reparameterization")
  rval
}

## Final-state fit in the weighted DR basis.
calc_smooth_dr_fit <- function(dr, lambda, native = TRUE)
{
  if(isTRUE(native) && is.double(dr$T))
    return(.Call(C_calc_smooth_dr_fit, dr, as.numeric(lambda), PACKAGE = "gamlss2"))
  q <- 1 / (1 + as.numeric(lambda[1L]) * dr$d)
  list(coefficients = drop(dr$T %*% (dr$g * q)), edf = sum(dr$h * q),
    vcov = tcrossprod(sweep(dr$T, 2L, sqrt(q), "*")))
}

## Quadratic smoothing criterion and its partial derivatives in RSS and EDF.
smooth.construct_score <- function(rss, edf, n, K, criterion)
{
  switch(criterion,
    "gcv" = {
      den <- n - edf
      list(value = rss * n / den^2,
        rss = n / den^2, edf = 2 * rss * n / den^3)
    },
    "aic" = list(value = rss + 2 * edf, rss = 1, edf = 2),
    "gaic" = list(value = rss + K * edf, rss = 1, edf = K),
    "aicc" = {
      den <- n - edf - 1
      list(value = rss + 2 * edf + 2 * edf * (edf + 1) / den,
        rss = 1, edf = 2 +
          2 * ((2 * edf + 1) * den + edf * (edf + 1)) / den^2)
    },
    "bic" = list(value = rss + log(n) * edf, rss = 1, edf = log(n))
  )
}


## Scale initial smoothing parameters to the weighted crossproduct matrix.
smooth.construct_start <- function(XWX, penalties)
{
  p <- ncol(XWX)
  m <- length(penalties)
  if(!m)
    return(numeric())

  diagonal <- pmax(0, diag(XWX))
  penalized <- rep(FALSE, p)
  scaled.diagonal <- numeric(p)
  lambda <- rep(1.0, m)

  for(i in seq_len(m)) {
    S <- penalties[[i]]
    if(!is.matrix(S) || !identical(dim(S), c(p, p)))
      next
    absolute <- abs(S)
    largest <- max(absolute)
    if(!is.finite(largest) || largest <= 0)
      next
    tolerance <- .Machine$double.eps^0.8 * largest
    active <- rowMeans(absolute) > tolerance &
      colMeans(absolute) > tolerance & diag(absolute) > tolerance
    if(!any(active))
      next
    size.XWX <- mean(diagonal[active])
    size.S <- mean(diag(S)[active])
    if(is.finite(size.XWX) && size.XWX > 0 &&
        is.finite(size.S) && size.S > 0)
      lambda[i] <- size.XWX / size.S
    penalized <- penalized | active
    scaled.diagonal <- scaled.diagonal + lambda[i] * diag(S)
  }

  active <- penalized & is.finite(diagonal) & diagonal > 0 &
    is.finite(scaled.diagonal) & scaled.diagonal > 0
  if(any(active)) {
    retained <- function()
      mean(diagonal[active] / (diagonal[active] + scaled.diagonal[active]))
    while(retained() > 0.4 && max(lambda) < 1e+10) {
      lambda <- 10 * lambda
      scaled.diagonal <- 10 * scaled.diagonal
    }
    while(retained() < 0.4 && min(lambda) > 1e-10) {
      lambda <- lambda / 10
      scaled.diagonal <- scaled.diagonal / 10
    }
  }

  pmin(1e+10, pmax(1e-10, lambda))
}

## Fitting function for mgcv smooth terms.
smooth.construct_wfit <- function(x, z, w, y, eta, j, family, control, transfer, iter)
{
  ## Fixed smooths can omit S entirely. Normalize with exact lookup: x$S
  ## would otherwise partially match S.scale and return a numeric vector.
  if(is.null(x[["S", exact = TRUE]]))
    x$S <- list()

  ## Number of observations.
  n <- length(z)

  if(control$binning) {
    rw <- numeric(length(x$binning$nodups))
    rz <- numeric(length(x$binning$nodups))
    calc_Xe(x$binning$sorted.index, w, z, rw, rz, x$binning$order)
  }

  ## Pre compute matrices.
  cache <- x$.rs_cache
  reuse <- is.environment(cache) && !is.null(cache$XWX) &&
    identical(control$binning, cache$binning) &&
    identical(x$X, cache$X, num.eq = FALSE, single.NA = FALSE) &&
    identical(w, cache$w, num.eq = FALSE, single.NA = FALSE)
  zWz <- NULL
  if(control$binning) {
    XWz <- crossprod(x$X, rz)
    reuse <- reuse && identical(x$binning, cache$bins) &&
      identical(x$sparse_index, cache$sparse_index)
    XWX <- if(reuse) cache$XWX else calc_XWX(x$X, 1/rw, x$sparse_index)
  } else {
    use.symmetric <- is.double(x$X) && ncol(x$X) >= 10L &&
      all(is.finite(w)) && all(w >= 0)
    if(use.symmetric && all(w == w[1L])) {
      ## Constant weights only rescale X'X. Keep its separate design key
      ## so changes in the weight level still reuse the unweighted matrix.
      reuse.XTX <- is.environment(cache) && !is.null(cache$XTX) &&
        identical(x$X, cache$XTX.X, num.eq = FALSE, single.NA = FALSE)
      reuse <- reuse && reuse.XTX
      crossproducts <- calc_XWXz(x$X, rep(1.0, n), z,
        if(reuse.XTX) cache$XTX else NULL)
      XWX <- if(reuse) cache$XWX else w[1L] * crossproducts$XWX
      XWz <- w[1L] * crossproducts$XWz
      zWz <- w[1L] * crossproducts$zWz
      if(is.environment(cache) && !reuse.XTX) {
        cache$XTX.X <- x$X
        cache$XTX <- crossproducts$XWX
      }
    } else if(use.symmetric) {
      ## The native kernel uses the same blocks and summation order for
      ## X'Wz and z'Wz with or without a cached X'WX.
      crossproducts <- calc_XWXz(x$X, w, z, if(reuse) cache$XWX else NULL)
      XWX <- crossproducts$XWX
      XWz <- crossproducts$XWz
      zWz <- crossproducts$zWz
    } else {
      ## Preserve the previous behavior for non-finite or negative weights.
      XW <- if(reuse) cache$XW else x$X * w
      XWX <- if(reuse) cache$XWX else crossprod(XW, x$X)
      XWz <- crossprod(XW, z)
    }
  }
  if(is.environment(cache) && !reuse) {
    cache$dr <- NULL
    cache$gradient.root <- NULL
    cache$X <- x$X
    cache$w <- w
    cache$binning <- control$binning
    cache$XWX <- XWX
    cache$XW <- if(!control$binning && !use.symmetric) XW else NULL
    cache$bins <- if(control$binning) x$binning else NULL
    cache$sparse_index <- if(control$binning) x$sparse_index else NULL
  }
  S <- diag(1e-05, ncol(x$X))

  if(!is.null(x$control)) {
    control[names(x$control)] <- x$control
    if(!is.null(control$method))
      control$criterion <- tolower(control$method)
  }
  if(is.null(control$criterion))
    control$criterion <- "aicc"

  control$criterion <- tolower(control$criterion)
  ncv <- identical(control$criterion, "ncv")
  ## NCV neighbourhoods are term-local metadata. A scalar specifies an
  ## ordered lag; a list supplies one deletion neighbourhood per observation.
  ncv.neigh <- NULL
  ncv.target <- NULL
  ncv.one <- TRUE
  ncv.lag <- NULL
  if(ncv) {
    ncv.config <- if(is.list(x$xt)) x$xt$ncv else NULL
    if(is.null(ncv.config))
      ncv.config <- control$ncv

    if(!is.null(ncv.config)) {
      if(is.numeric(ncv.config) && !is.list(ncv.config)) {
        if(length(ncv.config) != 1L || !is.finite(ncv.config) ||
            ncv.config < 0 || ncv.config != floor(ncv.config))
          stop("NCV lag must be one non-negative integer")
        lag <- as.integer(min(ncv.config, n))
        ncv.lag <- lag
        if(lag > 0L) {
          ncv.neigh <- lapply(seq_len(n), function(i)
            seq.int(max(1L, i - lag), min(n, i + lag)))
          ncv.target <- pmin(lag + 1L, seq_len(n))
          ncv.one <- FALSE
        }
      } else if(is.list(ncv.config)) {
        if(length(ncv.config) != n)
          stop("NCV neighbourhood list must have length ", n)
        ncv.neigh <- vector("list", n)
        ncv.target <- integer(n)
        for(i in seq_len(n)) {
          a <- ncv.config[[i]]
          if(!is.numeric(a) || !length(a) || any(!is.finite(a)) ||
              any(a != floor(a)) || any(a < 1) || any(a > n))
            stop("invalid NCV neighbourhood for observation ", i)
          a <- unique(as.integer(a))
          target <- match(i, a)
          if(is.na(target))
            stop("NCV neighbourhood for observation ", i,
              " does not contain its target")
          ncv.neigh[[i]] <- a
          ncv.target[i] <- target
        }
        ncv.one <- all(lengths(ncv.neigh) == 1L)
      } else {
        stop("NCV neighbourhoods must be a list or a non-negative lag")
      }
    }

    if(any(!is.finite(w)) || any(w < 0))
      stop("NCV requires finite non-negative weights")
    ncv.sqrt.w <- sqrt(w)
    ncv.Xw <- if(control$binning)
      x$X[x$binning$match.index, , drop = FALSE] * ncv.sqrt.w else
      x$X * ncv.sqrt.w
    ncv.zw <- ncv.sqrt.w * z
  }

  ## Extra penalty for selection.
  if(isTRUE(control$termselect)) {
    df <- ncol(x$X)
    bml <- drop(solve(XWX + diag(1e-08, df), XWz))

    pen <- function(b) {
      ## A zero starting fit needs a finite selection penalty.
      A <- 1 / rep(max(sqrt(sum(b^2)), .Machine$double.eps), df) *
        1 / rep(max(sqrt(sum(bml^2)), .Machine$double.eps), df)
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
    if(is.null(control$start)) {
      lambdas <- 1.0
      ## Small one-dimensional bases benefit from a scale-aware start in
      ## the full likelihood optimizer. Larger bases and multiple penalties
      ## are normally reached with warm starts from the outer RS iterations.
      if(isTRUE(control$logLik) && !isTRUE(x$localML) &&
          length(x$S) == 1L && isTRUE(x$dim == 1L) &&
          ncol(x$X) <= 30L && is.null(x$sp) && !isTRUE(x$fixed)) {
        scaled <- smooth.construct_start(XWX, x$S)
        if(length(scaled) == 1L && is.finite(scaled) && scaled >= 10)
          lambdas <- scaled
      }
    } else {
      lambdas <- control$start
    }
  }
  lambdas <- rep(lambdas, length.out = length(x$S))

  ## Penalty for AIC.
  K <- if(is.null(control$K)) 2 else control$K

  ## Local ML check.
  localML <- isTRUE(x$localML) && length(x$S) == 1L &&
    !isTRUE(x$fixed) && is.null(x$sp) && any(x$S[[1L]] != 0) &&
    !(length(x$rank) == 1L && x$rank == 0)
  if(!localML) {
    if(control$criterion == "ml")
      control$criterion <- "aicc"
  }

  ## Choose the direct or Demmler-Reinsch solver. The native quadratic
  ## search amortizes DR setup even for a short, warm-started update.
  dr.mode <- control$demmler.reinsch
  if(is.null(dr.mode))
    dr.mode <- "auto"
  if(is.character(dr.mode) && length(dr.mode) == 1L)
    dr.mode <- tolower(dr.mode)
  if(!(identical(dr.mode, "auto") || identical(dr.mode, TRUE) ||
      identical(dr.mode, FALSE)))
    stop("'demmler.reinsch' must be one of 'auto', TRUE or FALSE")

  dr.eligible <- length(x$S) == 1L &&
      is.null(x$sp) &&
      !isTRUE(x$fixed) &&
      (iter[1L] > 0L || localML)

  dr_setup <- function() {
    dr.ridge <- if(control$criterion == "ml" && localML) 0 else 1e-05
    ## The metric and penalty determine T, d, h and M; the response only
    ## determines g. Crossproduct reuse alone is insufficient if S changes.
    cached <- if(is.environment(cache)) cache$dr else NULL
    if(!is.null(cached) &&
        identical(cached$S, x$S[[1L]], num.eq = FALSE, single.NA = FALSE) &&
        identical(cached$ridge, dr.ridge) &&
        identical(cached$rank, x$rank) &&
        identical(cached$native, !identical(control$native.wfit, FALSE))) {
      rval <- cached$decomposition
      rval$g <- drop(crossprod(rval$T, XWz))
      return(rval)
    }
    rval <- try(smooth.construct_dr(
      XWX = XWX, XWz = XWz, S = x$S[[1L]], ridge = dr.ridge,
      rank = x$rank, native = !identical(control$native.wfit, FALSE)
    ), silent = TRUE)
    if(inherits(rval, "try-error"))
      return(NULL)
    if(is.environment(cache)) {
      decomposition <- rval
      decomposition$g <- NULL
      cache$dr <- list(S = x$S[[1L]], ridge = dr.ridge, rank = x$rank,
        native = !identical(control$native.wfit, FALSE),
        decomposition = decomposition)
    }
    rval
  }

  dr <- NULL
  dr.attempted <- FALSE
  if(dr.eligible && (identical(dr.mode, TRUE) ||
      (localML && !identical(dr.mode, FALSE)) ||
      (identical(dr.mode, "auto") && !isTRUE(control$logLik) &&
        !identical(control$native.wfit, FALSE) &&
        !identical(control$analytic.gradient, FALSE)))) {
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
    ## Constraints can remove part of the penalty null space: for example,
    ## a centered second-order P-spline has one unpenalized coefficient.
    null.dim <- x$null.space.dim
    if(is.null(null.dim)) {
      if(length(x$rank) == 1L) {
        null.dim <- ncol(XWX) - x$rank
      } else {
        values <- eigen(x$S[[1L]], symmetric = TRUE, only.values = TRUE)$values
        null.dim <- sum(abs(values) <=
          sqrt(.Machine$double.eps) * max(abs(values)))
      }
    }

    ## Evaluate b'Sb in penalty coordinates to avoid cancellation near
    ## the null space. The penalty basis is independent of X, w and z.
    penalty.root <- NULL
    if(is.null(dr)) {
      cached <- if(is.environment(cache)) cache$ml.penalty.root else NULL
      if(!is.null(cached) && identical(cached$S, x$S[[1L]],
          num.eq = FALSE, single.NA = FALSE) && identical(cached$null, null.dim)) {
        penalty.root <- cached$decomposition
      } else {
        penalty.root <- eigen(x$S[[1L]], symmetric = TRUE)
        tol <- sqrt(.Machine$double.eps) * max(1, abs(penalty.root$values))
        if(any(penalty.root$values < -tol))
          stop("penalty matrix is not positive semi-definite")
        penalty.root$values[penalty.root$values < 0] <- 0
        if(null.dim > 0) {
          null <- seq.int(ncol(XWX) - null.dim + 1L, ncol(XWX))
          if(any(abs(penalty.root$values[null]) > tol))
            stop("penalty null space does not match its rank")
          penalty.root$values[null] <- 0
        }
        if(is.environment(cache))
          cache$ml.penalty.root <- list(S = x$S[[1L]], null = null.dim,
            decomposition = penalty.root)
      }
    }

    if(!identical(control$native.wfit, FALSE) && is.double(x$X) &&
        is.double(XWX) && is.double(x$S[[1L]])) {
      rval <- try(.Call(C_calc_smooth_ml, x$X, as.numeric(z), as.numeric(w),
        XWX, as.numeric(XWz), x$S[[1L]], as.numeric(lambdas),
        as.numeric(null.dim), if(control$binning) x$binning$match.index else NULL,
        dr, penalty.root, PACKAGE = "gamlss2"), silent = TRUE)
      if(!inherits(rval, "try-error"))
        return(rval)
    }

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
      if(is.null(dr)) {
        alpha <- drop(crossprod(penalty.root$vectors, b))
        tau2 <- sum(penalty.root$values * alpha^2) / (edf - null.dim)
      } else {
        alpha <- dr$g / (1 + lambdas * dr$d)
        tau2 <- sum(dr$d * alpha^2) / (edf - null.dim)
      }

      if(tau2 < 1e-07) tau2 <- 1e-07
      lambdas.old <- lambdas
      lambdas <- sig2/tau2
      if(lambdas < 1e-07) lambdas <- 1e-07
      if(lambdas > 1e+07) lambdas <- 1e+07
      if(abs(lambdas - lambdas.old) < 1e-07 || lambdas > 1e+10) break
    }

    ## The last update changes lambda after b, EDF and P were computed.
    ## Return a complete fit at the reported lambda for either solver.
    if(is.null(dr)) {
      P <- try(chol2inv(chol(XWX + lambdas * x$S[[1L]])), silent = TRUE)
      if(inherits(P, "try-error"))
        P <- solve(XWX + lambdas * x$S[[1L]])
      b <- drop(P %*% XWz)
      edf <- sum(XWX * P)
    } else {
      drs <- calc_smooth_dr_fit(dr, lambdas, native = !identical(control$native.wfit, FALSE))
      b <- drs$coefficients
      P <- drs$vcov
      edf <- drs$edf
    }
    fit <- drop(x$X %*% b)
    if(control$binning)
      fit <- fit[x$binning$match.index]

    return(list("coefficients" = b, "fitted.values" = fit, "edf" = edf,
      "lambdas" = lambdas, "vcov" = P, "df" = n - edf))
  } else {
    if(is.null(zWz))
      zWz <- sum(w * z^2)
    criterion.evaluations <- 0L
    criterion.levels <- c("gcv", "aic", "gaic", "aicc", "bic")
    criterion.code <- match(tolower(control$criterion), criterion.levels)

    native.penalties <- all(vapply(x$S, function(Sk) {
      is.matrix(Sk) && is.double(Sk) &&
        identical(dim(Sk), dim(XWX))
    }, logical(1L)))
    use.native.wfit <-
      !identical(control$native.wfit, FALSE) &&
      !is.na(criterion.code) &&
      is.matrix(XWX) && is.double(XWX) &&
      native.penalties

    use.gradient <- !identical(control$analytic.gradient, FALSE) &&
      !isTRUE(control$logLik) && !is.na(criterion.code)
    quadratic.state <- NULL
    gradient.root <- if(is.environment(cache)) cache$gradient.root else NULL

    ## Share the objective and gradient state at each lambda. In DR
    ## coordinates B becomes I - ridge * T'T: retain that ridge correction
    ## in RSS as well as in EDF, so the objective is unchanged.
    quadratic_eval <- function(l, gradient = FALSE) {
      if(is.null(quadratic.state) ||
          !identical(l, quadratic.state$lambda) ||
          !identical(!is.null(dr), quadratic.state$dr)) {
        if(!is.null(dr) && use.native.wfit) {
          state <- .Call(C_calc_smooth_dr_eval, dr, as.numeric(l),
            as.numeric(zWz), as.numeric(n), as.numeric(K),
            as.integer(criterion.code), PACKAGE = "gamlss2")
        } else if(!is.null(dr)) {
          t <- as.numeric(l[1L]) * dr$d
          q <- 1 / (1 + t)
          if(any(!is.finite(q)) || any(q <= 0))
            stop("invalid smoothing parameter in Demmler-Reinsch fit")
          alpha <- dr$g * q
          Ma <- drop(dr$M %*% alpha)
          rss <- zWz - 2 * sum(dr$g * alpha) + sum(alpha^2) -
            dr$ridge * sum(alpha * Ma)
          edf <- sum(dr$h * q)
          alpha1 <- -alpha * t * q
          rss1 <- 2 * sum((alpha - dr$g - dr$ridge * Ma) * alpha1)
          edf1 <- -sum(dr$h * q^2 * t)
          score <- smooth.construct_score(rss, edf, n, K,
            criterion.levels[criterion.code])
          state <- list(value = score$value, rss = rss, edf = edf,
            gradient = score$rss * rss1 + score$edf * edf1)
        } else {
          if(use.native.wfit) {
            state <- calc_smooth_wfit(XWX, XWz, x$S, l, 1e-05,
              zWz, n, K, criterion.code,
              gradient = use.gradient && length(x$S) == 1L, state = TRUE)
          } else {
            A <- XWX + S
            for(k in seq_along(x$S))
              A <- A + l[k] * x$S[[k]]
            R <- chol(A)
            b <- drop(backsolve(R, forwardsolve(t(R), XWz)))
            P <- chol2inv(R)
            state <- list(coefficients = b, vcov = P, edf = sum(XWX * P))
            b <- state$coefficients
            state$residual <- drop(XWX %*% b) - drop(XWz)
            rss <- zWz - 2 * sum(b * XWz) + sum(b * (XWX %*% b))
            state$rss <- rss
            state$score <- smooth.construct_score(rss, state$edf, n, K,
              criterion.levels[criterion.code])
            state$value <- state$score$value
          }
        }
        ## Crossproduct subtraction loses precision for nearly exact fits.
        ## Calculate RSS and its score residuals from observations instead.
        if(zWz > 0 && state$rss <= sqrt(.Machine$double.eps) * zWz) {
          if(is.null(dr)) {
            stable <- calc_smooth_residual(x$X, state$coefficients, z, w,
              if(control$binning) x$binning$match.index else NULL,
              native = use.native.wfit)
            state$residual <- stable$residual
            state$gradient <- NULL
          } else {
            t <- as.numeric(l[1L]) * dr$d
            q <- 1 / (1 + t)
            alpha <- dr$g * q
            b <- drop(dr$T %*% alpha)
            b1 <- drop(dr$T %*% (-alpha * t * q))
            stable <- calc_smooth_residual(x$X, b, z, w,
              if(control$binning) x$binning$match.index else NULL, b1,
              native = use.native.wfit)
          }
          state$rss <- stable$rss
          state$score <- smooth.construct_score(state$rss, state$edf, n, K,
            criterion.levels[criterion.code])
          state$value <- state$score$value
          if(!is.null(dr))
            state$gradient <- state$score$rss * stable$rss1 +
              state$score$edf * (-sum(dr$h * q^2 * t))
        }
        state$lambda <- l
        state$dr <- !is.null(dr)
        quadratic.state <<- state
      }
      if(gradient && is.null(quadratic.state$gradient)) {
        state <- quadratic.state
        if(use.native.wfit) {
          if(is.null(gradient.root) && length(x$S) > 1L && ncol(XWX) >= 10L) {
            gradient.root <<- try(chol(XWX), silent = TRUE)
            if(is.environment(cache))
              cache$gradient.root <- gradient.root
          }
          root <- if(inherits(gradient.root, "try-error")) NULL else gradient.root
          state$gradient <- calc_smooth_wfit_gradient(XWX, XWz, x$S, l,
            state, 1e-05, n, K, criterion.code, root = root)
        } else {
          P <- state$vcov
          PBP <- P %*% XWX %*% P
          state$gradient <- vapply(seq_along(x$S), function(k) {
            b1 <- -l[k] * drop(P %*% (x$S[[k]] %*% state$coefficients))
            rss1 <- 2 * sum(b1 * state$residual)
            edf1 <- -l[k] * sum(PBP * x$S[[k]])
            state$score$rss * rss1 + state$score$edf * edf1
          }, numeric(1L))
        }
        quadratic.state <<- state
      }
      quadratic.state
    }

    ncv_score <- function(b, Q, gradient = FALSE, P = NULL) {
      e <- ncv.zw - drop(ncv.Xw %*% b)
      if(any(!is.finite(e)) || any(!is.finite(Q)))
        return(Inf)

      if(!is.null(ncv.lag) && ncv.lag > 0L && gradient)
        return(.Call(C_calc_ncv_lag_gradient, Q, e, ncv.Xw, b, P,
          x$S, ncv.lag, PACKAGE = "gamlss2"))
      if(!is.null(ncv.lag))
        return(.Call(C_calc_ncv_lag, Q, e, ncv.lag, 0L,
          PACKAGE = "gamlss2"))

      if(ncv.one) {
        den <- 1 - rowSums(Q^2)
        if(any(!is.finite(den)) || any(den <= sqrt(.Machine$double.eps)))
          return(Inf)
        e <- e / den
        return(if(all(is.finite(e))) sum(e^2) else Inf)
      }

      ## For penalized WLS, deleted residuals are exactly
      ## (I - H[a, a])^-1 e[a], with H = Q Q'.
      value <- 0
      for(i in seq_len(n)) {
        a <- ncv.neigh[[i]]
        Qa <- Q[a, , drop = FALSE]
        Ra <- try(chol(diag(length(a)) - tcrossprod(Qa)), silent = TRUE)
        if(inherits(Ra, "try-error") || any(!is.finite(Ra)) ||
            any(diag(Ra)^2 <= sqrt(.Machine$double.eps)))
          return(Inf)
        u <- backsolve(Ra, forwardsolve(t(Ra), e[a]))
        target <- u[ncv.target[i]]
        if(!is.finite(target))
          return(Inf)
        value <- value + target^2
      }
      if(is.finite(value)) value else Inf
    }

    ncv_eval <- function(l, gradient = FALSE) {
      if(!is.null(dr)) {
        q <- 1 / (1 + as.numeric(l[1L]) * dr$d)
        if(any(!is.finite(q)) || any(q <= 0))
          return(Inf)
        alpha <- dr$g * q
        b <- drop(dr$T %*% alpha)
        Q <- ncv.Xw %*% sweep(dr$T, 2L, sqrt(q), "*")
        P <- dr$T %*% (q * t(dr$T))
      } else {
        Sl <- S
        for(k in seq_along(x$S))
          Sl <- Sl + l[k] * x$S[[k]]
        R <- try(chol(XWX + Sl), silent = TRUE)
        if(inherits(R, "try-error"))
          return(Inf)
        b <- drop(backsolve(R, forwardsolve(t(R), XWz)))
        Q <- t(forwardsolve(t(R), t(ncv.Xw)))
        P <- chol2inv(R)
      }
      ncv_score(b, Q, gradient = gradient, P = P)
    }

    ## Function to search for smoothing parameters using GCV etc.
    fl <- function(l, rf = FALSE) {
      if(!rf)
        criterion.evaluations <<- criterion.evaluations + 1L
      if(ncv && !rf)
        return(ncv_eval(l))
      if(!isTRUE(control$logLik) && !is.na(criterion.code) && !rf)
        return(quadratic_eval(l)$value)
      if(rf && !is.null(dr)) {
        drs <- calc_smooth_dr_fit(dr, l, native = use.native.wfit)
        b <- drs$coefficients
        edf <- drs$edf
        P <- drs$vcov
      } else if(rf && !is.null(quadratic.state) && !quadratic.state$dr &&
          identical(l, quadratic.state$lambda) &&
          !is.null(quadratic.state$vcov)) {
        b <- quadratic.state$coefficients
        edf <- quadratic.state$edf
        P <- quadratic.state$vcov
      } else if(!is.null(dr) && !rf) {
        drs <- dr_fit(l)
        edf <- drs$edf
        b <- drs$coefficients
      } else if(use.native.wfit) {
        native.fit <- calc_smooth_wfit(
          XWX = XWX,
          XWz = XWz,
          penalties = x$S,
          lambda = l,
          ridge = 1e-05,
          zWz = zWz,
          n = n,
          K = K,
          criterion = criterion.code,
          final = rf || isTRUE(control$logLik)
        )
        if(!rf && !isTRUE(control$logLik))
          return(native.fit)

        b <- native.fit$coefficients
        edf <- native.fit$edf
        P <- native.fit$vcov
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

    ## Check for fx = TRUE. Fully fixed tensor smooths may not have a scalar
    ## fixed flag, but without penalties there is nothing to optimize.
    if(isTRUE(x$fixed) || !length(x$S)) {
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
          gradient = if(use.gradient && !ncv) function(rho)
            quadratic_eval(exp(rho), gradient = TRUE)$gradient else
            if(ncv && !is.null(ncv.lag) && length(x$S)) function(rho) {
              out <- ncv_eval(exp(rho), gradient = TRUE)
              as.numeric(out[[2L]]) * exp(rho)
            } else NULL,
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

  ## Adaptive term selection changes the penalty during fitting. Retain only
  ## this otherwise unrecoverable final effective penalty.
  if(isTRUE(control$termselect)) {
    rval$penalty <- matrix(0, ncol(x$X), ncol(x$X))
    for(k in seq_along(x$S))
      rval$penalty <- rval$penalty +
        rval$lambdas[k] * x$S[[k]]
  }

  rval$transfer <- list(
    "lambdas" = rval$lambdas,
    "coefficients" = rval$coefficients,
    "criterion.evaluations" = criterion.evaluations,
    "demmler.reinsch" = !is.null(dr),
    "names" = colnames(x$X)
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
