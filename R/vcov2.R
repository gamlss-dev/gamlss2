## Joint covariance matrix for gamlss2 models -----------------------------
##
## This implementation is self-contained and does not call the superseded
## covariance helpers in inference.R.
##
## method = "joint" builds the complete observed likelihood information by
## differentiating the linked-scale family scores. Legacy family Hessian
## callbacks can be Fisher/working approximations and therefore are not, in
## general, the observed Hessian (in particular for nonlinear links).
##
## method = "working" uses final parameter-wise working curvature and zero
## cross-parameter blocks.
##
## method = "numeric" computes an independent coefficient-dimensional Hessian
## of the unpenalized negative log-likelihood and then adds the same penalties.
##
## full = FALSE is applied only after the complete joint inverse is formed.
##
## For special terms, a stored coefficient-linear design is used when
## available. Compact specials can be reconstructed through special_predict()
## only after deterministic checks verify affine dependence on coefficients.

vcov.gamlss2 <- function(object,
  type = c("vcov", "cor", "se", "coef"), full = FALSE,
  method = c("joint", "working", "numeric"), ...)
{
  type <- match.arg(type)
  method <- match.arg(method)

  ## utilities -----------------------------------------------------------

  .stop <- function(...)
    stop(..., call. = FALSE)

  .finite <- function(x)
    is.numeric(x) && all(is.finite(x))

  .matrix <- function(x)
  {
    if(inherits(x, "Matrix"))
      x <- as.matrix(x)
    if(!is.matrix(x) || !is.numeric(x))
      return(NULL)
    return(x)
  }

  .same <- function(x, y, tol = 1e-7)
  {
    if(length(x) != length(y))
      return(FALSE)
    ii <- is.na(x) & is.na(y)
    if(any(xor(is.na(x), is.na(y))))
      return(FALSE)
    if(all(ii))
      return(TRUE)
    x <- x[!ii]
    y <- y[!ii]
    sc <- max(1, abs(x), abs(y))
    return(max(abs(x - y)) <= tol * sc)
  }

  .fixed <- function(name, pos = NULL)
  {
    z <- object$control$fixed
    if(is.null(z) || !length(z))
      return(FALSE)
    if(!is.null(names(z)) && name %in% names(z))
      z <- z[[name]]
    else if(!is.null(pos) && length(z) >= pos)
      z <- z[[pos]]
    else
      return(FALSE)
    return(length(z) == 1L && !is.na(z) && as.logical(z))
  }

  .get_named <- function(x, name, pos = NULL)
  {
    if(is.null(x))
      return(NULL)
    if(!is.null(names(x)) && name %in% names(x))
      return(x[[name]])
    if(!is.null(pos) && length(x) >= pos)
      return(x[[pos]])
    return(NULL)
  }

  .coef_vector <- function(x)
  {
    if(is.null(x))
      return(numeric())
    if(is.numeric(x) && is.null(dim(x)))
      return(x)
    if(is.matrix(x) && (ncol(x) == 1L || nrow(x) == 1L))
      return(as.numeric(x))
    if(is.list(x)) {
      z <- unlist(x, recursive = TRUE, use.names = TRUE)
      if(is.numeric(z))
        return(z)
    }
    return(NULL)
  }

  .eta_list <- function(x, names)
  {
    out <- vector("list", length(names))
    names(out) <- names

    if(is.data.frame(x) || is.matrix(x)) {
      if(nrow(x) < 1L)
        .stop("empty fitted predictors in 'object$fitted.values'.")
      for(j in seq_along(names)) {
        if(!is.null(colnames(x)) && names[j] %in% colnames(x))
          out[[j]] <- as.numeric(x[, names[j]])
        else if(ncol(x) >= j)
          out[[j]] <- as.numeric(x[, j])
        else
          .stop("cannot match fitted predictors to family parameters.")
      }
      return(out)
    }

    if(is.list(x)) {
      for(j in seq_along(names)) {
        z <- .get_named(x, names[j], j)
        if(is.null(z))
          .stop("cannot match fitted predictors to family parameters.")
        out[[j]] <- as.numeric(z)
      }
      return(out)
    }

    if(length(names) == 1L && is.numeric(x)) {
      out[[1L]] <- as.numeric(x)
      return(out)
    }

    .stop("unsupported 'object$fitted.values' representation.")
  }

  .response <- function(n)
  {
    y <- object$y

    if(is.null(y) && !is.null(object$model))
      y <- tryCatch(stats::model.response(object$model),
        error = function(e) NULL)

    if(is.null(y)) {
      mf <- tryCatch(stats::model.frame(object), error = function(e) NULL)
      if(!is.null(mf))
        y <- tryCatch(stats::model.response(mf), error = function(e) NULL)
    }

    if(is.null(y))
      .stop(paste0(
        "the response is not available. Models fitted with 'light = TRUE' ",
        "must be refitted before joint covariance can be computed."
      ))

    ny <- if(is.matrix(y) || is.data.frame(y)) nrow(y) else length(y)
    if(ny != n)
      .stop("response and fitted predictors have incompatible lengths.")

    return(y)
  }

  .prior_weights <- function(n)
  {
    w <- object$weights

    if(is.null(w) && !is.null(object$model))
      w <- tryCatch(stats::model.weights(object$model),
        error = function(e) NULL)

    if(is.null(w))
      w <- rep(1, n)

    w <- as.numeric(w)

    if(length(w) != n || any(!is.finite(w)) || any(w < 0))
      .stop("invalid prior observation weights in fitted model.")

    return(w)
  }

  .linear_matrix <- function(name, b, n)
  {
    p <- length(b)
    if(!p)
      return(matrix(0, n, 0L))

    bn <- names(b)
    b0 <- as.numeric(b)
    b0[!is.finite(b0)] <- 0

    fl <- .get_named(object$fitted.linear, name)
    fv <- if(is.list(fl)) fl$fitted.values else NULL
    if(!is.null(fv))
      fv <- as.numeric(fv)

    .align <- function(X)
    {
      X <- .matrix(X)
      if(is.null(X) || nrow(X) != n)
        return(NULL)

      if(ncol(X) == p) {
        if(!is.null(bn) && !is.null(colnames(X))) {
          i <- match(bn, colnames(X))
          if(!anyNA(i))
            X <- X[, i, drop = FALSE]
          else {
            pbn <- paste0(name, ".", bn)
            i <- match(pbn, colnames(X))
            if(!anyNA(i))
              X <- X[, i, drop = FALSE]
          }
        }
      } else {
        if(is.null(bn) || is.null(colnames(X)))
          return(NULL)

        i <- match(bn, colnames(X))
        if(anyNA(i)) {
          pbn <- paste0(name, ".", bn)
          i <- match(pbn, colnames(X))
        }
        if(anyNA(i))
          return(NULL)
        X <- X[, i, drop = FALSE]
      }

      if(ncol(X) != p || any(!is.finite(X)))
        return(NULL)

      ## protect against changed contrasts or an incorrectly selected matrix
      if(!is.null(fv) && length(fv) == n) {
        z <- drop(X %*% b0)
        sc <- max(1, abs(z), abs(fv))
        if(max(abs(z - fv)) > 1e-6 * sc)
          return(NULL)
      }

      return(X)
    }

    ## the linear terms object is the most reliable representation
    xt <- .get_named(object$xterms, name)
    if(!is.null(xt) && !is.null(object$model)) {
      X <- tryCatch(
        stats::model.matrix(xt, data = object$model),
        error = function(e) NULL
      )
      X <- .align(X)
      if(!is.null(X))
        return(X)
    }

    ## use stored parameter-specific linear fit information if available
    if(is.list(fl)) {
      for(nm in c("X", "x", "model.matrix")) {
        X <- .align(fl[[nm]])
        if(!is.null(X))
          return(X)
      }
    }

    ## some optimizers store one matrix for each distribution parameter
    if(is.list(object$x)) {
      X <- .align(.get_named(object$x, name))
      if(!is.null(X))
        return(X)
    }

    ## finally try the public model.matrix() method. This is only
    ## accepted if its dimensions, names and fitted linear contribution can
    ## be verified against the stored fit.
    mm <- tryCatch(
      getS3method("model.matrix", "gamlss2", optional = TRUE),
      error = function(e) NULL
    )

    if(is.function(mm)) {
      aa <- list(
        list(object = object, model = name),
        list(object = object, parameter = name),
        list(object = object)
      )

      for(a in aa) {
        Z <- tryCatch(do.call(mm, a), error = function(e) NULL)

        if(is.list(Z) && !is.data.frame(Z)) {
          zz <- .get_named(Z, name)
          if(!is.null(zz))
            Z <- zz
        }

        X <- .align(Z)
        if(!is.null(X))
          return(X)
      }
    }

    ## stored full model matrix
    X <- .align(object$x)
    if(!is.null(X))
      return(X)

    .stop(paste0(
      "cannot reconstruct and verify the linear model matrix for parameter '",
      name, "'."
    ))
  }

  .special_source <- function(name, term)
  {
    z <- object$specials
    if(is.null(z))
      return(NULL)

    if(is.list(z) && !is.null(names(z)) && term %in% names(z))
      return(z[[term]])

    if(is.list(z) && !is.null(names(z)) && name %in% names(z)) {
      zz <- z[[name]]
      if(is.list(zz) && !is.null(names(zz)) && term %in% names(zz))
        return(zz[[term]])
    }

    if(is.list(z)) {
      hit <- which(vapply(z, function(a) {
        is.list(a) && !is.null(a$label) &&
          identical(as.character(a$label), term)
      }, logical(1L)))
      if(length(hit) == 1L)
        return(z[[hit]])
    }

    return(NULL)
  }

  .special_matrix <- function(name, term, fit, source, b, n)
  {
    q <- length(b)
    if(!q)
      .stop(paste0(
        "special term '", term, "' in parameter '", name,
        "' has no coefficient vector."
      ))

    candidates <- list()
    for(z in list(fit, source)) {
      if(is.list(z)) {
        for(nm in c("X", "x", "model.matrix", "design"))
          if(!is.null(z[[nm]]))
            candidates[[length(candidates) + 1L]] <- z[[nm]]

        for(nn in c("smooth", "model", "fit")) {
          zz <- z[[nn]]
          if(is.list(zz)) {
            for(nm in c("X", "x", "model.matrix", "design"))
              if(!is.null(zz[[nm]]))
                candidates[[length(candidates) + 1L]] <- zz[[nm]]
          }
        }
      }
    }

    fv <- fit$fitted.values
    if(!is.null(fv))
      fv <- as.numeric(fv)

    b0 <- as.numeric(b)
    b0[!is.finite(b0)] <- 0

    .ok <- function(X)
    {
      X <- .matrix(X)
      if(!is.null(X) && nrow(X) != n) {
        binning <- if(is.list(source)) source$binning else NULL
        mi <- if(is.list(binning)) binning$match.index else NULL
        if(!is.null(mi) && length(mi) == n &&
          all(is.finite(mi)) && all(mi >= 1L) && all(mi <= nrow(X)))
          X <- X[as.integer(mi), , drop = FALSE]
      }
      if(is.null(X) || nrow(X) != n || ncol(X) != q ||
        any(!is.finite(X)))
        return(NULL)

      if(is.null(fv) || length(fv) != n)
        return(X)

      z <- drop(X %*% b0)
      sc <- max(1, abs(z), abs(fv))
      if(max(abs(z - fv)) <= 1e-6 * sc)
        return(X)

      ## smooth terms are centered before entering the additive predictor.
      Xc <- sweep(X, 2L, colMeans(X), "-")
      z <- drop(Xc %*% b0)
      sc <- max(1, abs(z), abs(fv))
      if(max(abs(z - fv)) <= 1e-6 * sc)
        return(Xc)

      return(NULL)
    }

    for(X in candidates) {
      X <- .ok(X)
      if(!is.null(X))
        return(X)
    }

    ## Some coefficient-linear specials store a compact representation
    ## instead of the complete training design matrix. Reconstruct the
    ## coefficient Jacobian through special_predict() and accept it only if
    ## deterministic perturbation checks confirm coefficient-linearity.
    data <- object$model
    spredict <- tryCatch(
      getFromNamespace("special_predict", "gamlss2"),
      error = function(e) NULL
    )

    if(!is.null(data) && is.function(spredict)) {
      .predict <- function(cf)
      {
        z <- fit
        z$coefficients <- cf

        v <- tryCatch(
          spredict(z, data = data, se.fit = FALSE),
          error = function(e) NULL
        )

        if(is.data.frame(v)) {
          if("fit" %in% names(v))
            v <- v$fit
          else if(ncol(v) == 1L)
            v <- v[[1L]]
        }

        if(is.matrix(v) && ncol(v) == 1L)
          v <- v[, 1L]

        v <- as.numeric(v)

        if(length(v) != n || any(!is.finite(v)))
          return(NULL)

        return(v)
      }

      f0 <- .predict(b0)

      if(!is.null(f0)) {
        X <- matrix(NA_real_, n, q)

        for(j in seq_len(q)) {
          h <- .Machine$double.eps^(1 / 3) * (1 + abs(b0[j]))
          bp <- bm <- b0
          bp[j] <- bp[j] + h
          bm[j] <- bm[j] - h

          fp <- .predict(bp)
          fm <- .predict(bm)

          if(is.null(fp) || is.null(fm)) {
            X <- NULL
            break
          }

          X[, j] <- (fp - fm) / (2 * h)
        }

        if(!is.null(X) && all(is.finite(X))) {
          ## A local Jacobian is not sufficient for nonlinear specials.
          ## Verify affine coefficient dependence in several deterministic
          ## directions before using X as a coefficient design matrix.
          ok <- TRUE
          direction <- list(
            rep(1, q),
            rep(c(-1, 1), length.out = q),
            seq_len(q) / max(1, q)
          )

          for(dd in direction) {
            dd <- as.numeric(dd)
            scd <- max(1, sqrt(sum(dd^2)))
            dd <- dd / scd
            h <- 1e-4 * (1 + max(abs(b0)))
            f1 <- .predict(b0 + h * dd)

            if(is.null(f1)) {
              ok <- FALSE
              break
            }

            target <- f0 + h * drop(X %*% dd)
            sc <- max(1, abs(f0), abs(f1), abs(target))

            if(max(abs(f1 - target)) > 1e-6 * sc) {
              ok <- FALSE
              break
            }
          }

          if(ok) {
            X <- .ok(X)
            if(!is.null(X))
              return(X)
          }
        }
      }
    }

    .stop(paste0(
      "cannot obtain a verified coefficient-linear design matrix for ",
      "special term '", term, "' in parameter '", name,
      "'. Nonlinear specials are not supported by joint covariance."
    ))
  }

  .special_penalty <- function(name, term, fit, source, q)
  {
    ## Active inequality constraints change the local coefficient geometry.
    ## A plain inverse of H + P is not the covariance conditional on an
    ## active constraint set. Do not silently treat such terms as ordinary
    ## quadratic smooths.
    cl <- class(fit)
    if(any(grepl("^ms[.]|^ms$", cl)) ||
      !is.null(fit$active) || !is.null(fit$active.set))
      .stop(paste0(
        "special term '", term, "' in parameter '", name,
        "' uses active constraints. A plain H + P covariance is not ",
        "valid for this term."
      ))

    ## Prefer an explicitly stored final penalty if the fitted term provides
    ## one. It must already be on the fitted coefficient scale.
    for(nm in c("penalty", "P.final", "penalty.matrix")) {
      z <- if(is.list(fit)) fit[[nm]] else NULL
      z <- .matrix(z)

      if(!is.null(z)) {
        if(any(dim(z) != q) || any(!is.finite(z)))
          .stop(paste0(
            "invalid fitted penalty matrix for special term '", term, "'."
          ))

        z <- 0.5 * (z + t(z))
        ev <- eigen(z, symmetric = TRUE, only.values = TRUE)$values
        cut <- 1e-9 * max(1, max(abs(ev)))

        if(any(ev < -cut))
          .stop(paste0(
            "fitted penalty matrix for special term '", term,
            "' is not positive semidefinite."
          ))

        return(z)
      }
    }

    ## Otherwise reconstruct the final quadratic penalty from the unscaled
    ## penalty matrices and the fitted smoothing parameter(s).
    S <- NULL
    if(is.list(source) && !is.null(source$S))
      S <- source$S
    if(is.null(S) && is.list(fit) && !is.null(fit$S))
      S <- fit$S

    if(is.null(S) && is.list(source) && isTRUE(source$fixed))
      return(matrix(0, q, q))

    if(is.null(S)) {
      lambda0 <- fit$lambdas
      if(is.null(lambda0))
        lambda0 <- fit$lambda
      if(!is.null(lambda0)) {
        lambda0 <- as.numeric(lambda0)
        if(length(lambda0) && all(is.finite(lambda0)) && all(lambda0 == 0))
          return(matrix(0, q, q))
      }
      .stop(paste0(
        "special term '", term, "' in parameter '", name,
        "' has coefficients but no verifiable quadratic penalty ",
        "representation."
      ))
    }

    ## Fixed-basis smooths legitimately have no penalty matrices.
    if(!length(S))
      return(matrix(0, q, q))

    if(is.matrix(S) || inherits(S, "Matrix"))
      S <- list(S)

    if(!is.list(S) || !length(S))
      .stop(paste0(
        "invalid penalty representation for special term '", term, "'."
      ))

    for(j in seq_along(S)) {
      S[[j]] <- .matrix(S[[j]])

      if(is.null(S[[j]]) || any(dim(S[[j]]) != q) ||
        any(!is.finite(S[[j]])))
        .stop(paste0(
          "invalid penalty matrix for special term '", term, "'."
        ))

      S[[j]] <- 0.5 * (S[[j]] + t(S[[j]]))
      ev <- eigen(S[[j]], symmetric = TRUE, only.values = TRUE)$values
      cut <- 1e-9 * max(1, max(abs(ev)))

      if(any(ev < -cut))
        .stop(paste0(
          "penalty matrix for special term '", term,
          "' is not positive semidefinite."
        ))
    }

    lambda <- fit$lambdas
    if(is.null(lambda))
      lambda <- fit$lambda

    ## For a genuinely fixed smooth, the constructor value can be used.
    ## Negative mgcv-style values denote an unknown smoothing parameter and
    ## must not be interpreted as fitted lambda.
    if(is.null(lambda) && is.list(source) && !is.null(source$sp)) {
      z <- as.numeric(source$sp)
      if(length(z) == length(S) && all(is.finite(z)) && all(z >= 0))
        lambda <- z
    }

    if(is.null(lambda)) {
      .stop(paste0(
        "cannot find the fitted smoothing parameter(s) for special term '",
        term, "' in parameter '", name, "'."
      ))
    }

    lambda <- as.numeric(lambda)

    if(length(lambda) != length(S) ||
      any(!is.finite(lambda)) || any(lambda < 0))
      .stop(paste0(
        "smoothing parameter(s) and penalty matrices do not match for ",
        "special term '", term, "'."
      ))

    P <- matrix(0, q, q)
    for(j in seq_along(S))
      P <- P + lambda[j] * S[[j]]

    return(P)
  }

  .parameter_block <- function(name, pos, n)
  {
    b <- .get_named(object$coefficients, name, pos)
    if(is.null(b))
      b <- numeric()
    b <- as.numeric(b)
    names(b) <- names(.get_named(object$coefficients, name, pos))

    X <- .linear_matrix(name, b, n)
    fl <- .get_named(object$fitted.linear, name, pos)
    p0 <- length(b)

    cn <- if(p0) {
      z <- names(b)
      if(is.null(z) || any(!nzchar(z)))
        z <- paste0("linear", seq_len(p0))
      paste0(name, ".", z)
    } else character()

    P <- matrix(0, p0, p0)
    if(isTRUE(object$control$ridge) && p0) {
      lambda <- if(is.list(fl)) fl$penalty else NULL
      if(is.null(lambda))
        .stop(paste0(
          "final linear ridge penalty is unavailable for parameter '",
          name, "'."
        ))
      lambda <- as.numeric(lambda)
      if(length(lambda) != 1L || !is.finite(lambda) || lambda < 0)
        .stop(paste0(
          "invalid final linear ridge penalty for parameter '", name, "'."
        ))
      penalized <- if(is.null(names(b))) rep(TRUE, p0) else
        names(b) != "(Intercept)"
      P <- diag(lambda * penalized, p0)
    }
    linear <- seq_len(p0)
    term <- rep("linear", p0)

    fs <- .get_named(object$fitted.specials, name, pos)

    if(!is.null(fs) && length(fs)) {
      if(is.null(names(fs)))
        .stop(paste0(
          "fitted special terms for parameter '", name,
          "' are not named."
        ))

      for(tt in names(fs)) {
        fit <- fs[[tt]]

        if(!is.list(fit) || is.null(fit$coefficients))
          .stop(paste0(
            "special term '", tt, "' in parameter '", name,
            "' has no coefficient-linear representation. Joint covariance ",
            "does not support this nonlinear special."
          ))

        bs <- .coef_vector(fit$coefficients)
        if(is.null(bs) || !length(bs))
          .stop(paste0(
            "invalid coefficient vector for special term '", tt, "'."
          ))

        source <- .special_source(name, tt)
        Xs <- .special_matrix(name, tt, fit, source, bs, n)
        Ps <- .special_penalty(name, tt, fit, source, length(bs))

        old <- ncol(X)
        X <- cbind(X, Xs)
        b <- c(b, bs)

        sn <- names(bs)
        if(is.null(sn) || any(!nzchar(sn)))
          sn <- seq_along(bs)
        cn <- c(cn, paste0(name, ".", tt, ".", sn))
        term <- c(term, rep(tt, length(bs)))

        Q <- matrix(0, nrow(P) + nrow(Ps), ncol(P) + ncol(Ps))
        if(length(P))
          Q[seq_len(nrow(P)), seq_len(ncol(P))] <- P
        ii <- old + seq_len(nrow(Ps))
        Q[ii, ii] <- Ps
        P <- Q
      }
    }

    if(!length(b))
      return(NULL)

    if(nrow(X) != n || ncol(X) != length(b))
      .stop(paste0(
        "coefficient/design dimension mismatch for parameter '", name, "'."
      ))

    colnames(X) <- cn
    rownames(P) <- colnames(P) <- cn

    return(list(
      name = name,
      coefficients = b,
      X = X,
      P = P,
      linear = linear,
      term = term,
      names = cn,
      fixed = rep(.fixed(name, pos), length(b))
    ))
  }

  .canonical_coef <- function(beta, full, linear)
  {
    z <- tryCatch(
      stats::coef(object, full = full, dropall = FALSE),
      error = function(e) NULL
    )
    z <- .coef_vector(z)

    target <- if(full) beta else beta[linear]

    if(!is.null(z) && length(z) == length(target) &&
      .same(as.numeric(z), as.numeric(target))) {
      if(!is.null(names(z)) && all(nzchar(names(z))))
        names(target) <- names(z)
    }

    return(target)
  }

  .eval_hessian <- function(family, par, y, a, b, n)
  {
    nm1 <- if(identical(a, b)) a else paste(a, b, sep = ":")
    nm2 <- if(identical(a, b)) a else paste(b, a, sep = ":")

    f1 <- family$hessian[[nm1]]
    f2 <- if(identical(nm1, nm2)) NULL else family$hessian[[nm2]]

    if(!is.function(f1) && !is.function(f2))
      .stop(paste0(
        "family Hessian component '", nm1,
        "' is not available. The fitted family must be completed before ",
        "joint covariance is computed."
      ))

    .eval <- function(f, nm)
    {
      z <- f(par = par, y = y)

      if(length(z) == 1L)
        z <- rep(z, n)

      z <- as.numeric(z)

      if(length(z) != n || any(!is.finite(z)))
        .stop(paste0(
          "family Hessian component '", nm,
          "' returned invalid values."
        ))

      return(z)
    }

    if(is.function(f1))
      z1 <- .eval(f1, nm1)
    else
      z1 <- NULL

    if(is.function(f2))
      z2 <- .eval(f2, nm2)
    else
      z2 <- NULL

    ## If complete_family() exposes both orientations, they must represent
    ## the same mixed derivative. Check this rather than choosing silently.
    if(!is.null(z1) && !is.null(z2)) {
      sc <- max(1, abs(z1), abs(z2))
      if(max(abs(z1 - z2)) > 1e-7 * sc)
        .stop(paste0(
          "family cross Hessians '", nm1, "' and '", nm2,
          "' disagree at the fitted values."
        ))
    }

    if(!is.null(z1))
      return(z1)

    return(z2)
  }

  .observed_curvature <- function(family, eta, y, a, b, n)
  {
    if(!is.function(family$map2par))
      .stop("the fitted family has no predictor-to-parameter map.")

    .differentiate <- function(score.name, predictor.name)
    {
      score <- family$score[[score.name]]
      if(!is.function(score))
        .stop(paste0(
          "family score component '", score.name, "' is not available."
        ))

      step <- .Machine$double.eps^(1 / 3) *
        pmax(1, abs(eta[[predictor.name]]))
      upper <- lower <- eta
      upper[[predictor.name]] <- eta[[predictor.name]] + step
      lower[[predictor.name]] <- eta[[predictor.name]] - step

      su <- score(par = family$map2par(upper), y = y)
      sl <- score(par = family$map2par(lower), y = y)
      z <- -(as.numeric(su) - as.numeric(sl)) / (2 * step)

      if(length(z) != n || any(!is.finite(z)))
        .stop(paste0(
          "cannot obtain finite observed curvature for parameters '",
          score.name, "' and '", predictor.name, "'."
        ))

      z
    }

    z1 <- .differentiate(a, b)
    if(identical(a, b))
      return(z1)

    ## Both orientations are derivatives of the same scalar log likelihood.
    ## Averaging removes small finite-difference asymmetry, while a generous
    ## check catches incompatible custom family score functions.
    z2 <- .differentiate(b, a)
    sc <- pmax(1, abs(z1), abs(z2))
    err <- max(abs(z1 - z2) / sc)
    if(err > 1e-3)
      .stop(paste0(
        "linked family scores give inconsistent mixed curvature for '",
        a, "' and '", b, "' (maximum relative discrepancy ",
        format(err, scientific = TRUE), ")."
      ))

    0.5 * (z1 + z2)
  }

  .working_weights <- function(family, par, y, eta, name, n)
  {
    f <- family$update
    if(is.list(f) && !is.null(f[[name]]))
      f <- f[[name]]

    if(is.function(f)) {
      z <- tryCatch(
        f(par = par, y = y, eta = eta, which = name),
        error = function(e) NULL
      )

      if(is.list(z) && !is.null(z$weights)) {
        w <- as.numeric(z$weights)
        if(length(w) == 1L)
          w <- rep(w, n)
        if(length(w) == n && all(is.finite(w)) && all(w >= 0))
          return(w)
      }
    }

    ## completed families define the same working curvature through the
    ## diagonal Hessian when no optimized update method is available
    w <- .eval_hessian(family, par, y, name, name, n)

    if(any(w < 0))
      .stop(paste0(
        "negative working curvature for parameter '", name,
        "'. No positive family working weights are available."
      ))

    return(w)
  }

  .block_diag <- function(blocks, what = c("X", "P"))
  {
    what <- match.arg(what)

    if(what == "X")
      .stop("internal misuse of '.block_diag'.")

    p <- sum(vapply(blocks, function(z) nrow(z$P), integer(1L)))
    out <- matrix(0, p, p)
    k <- 0L
    for(z in blocks) {
      q <- nrow(z$P)
      ii <- k + seq_len(q)
      out[ii, ii] <- z$P
      k <- k + q
    }
    return(out)
  }

  .invert <- function(K, names, forced.na, fixed, noise = 0)
  {
    K <- .matrix(K)

    if(is.null(K) || nrow(K) != ncol(K) || any(!is.finite(K)))
      .stop("invalid joint penalized information matrix.")

    sc <- max(1, max(abs(K)))
    asym <- max(abs(K - t(K))) / sc
    if(asym > 1e-7)
      .stop(paste0(
        "joint penalized information is not symmetric; relative asymmetry = ",
        format(asym, scientific = TRUE), "."
      ))

    K <- 0.5 * (K + t(K))
    p <- nrow(K)
    if(length(forced.na) != p || length(fixed) != p)
      .stop("invalid coefficient status map.")

    ## Fixed coefficients have exactly zero covariance. Coefficients already
    ## marked as aliased by the fit retain NA rows and columns.
    V <- matrix(0, p, p)
    aliased <- forced.na & !fixed
    if(any(aliased))
      V[aliased, ] <- V[, aliased] <- NA_real_
    active <- which(!forced.na & !fixed)

    if(!length(active)) {
      dimnames(V) <- list(names, names)
      return(V)
    }

    Ka <- K[active, active, drop = FALSE]
    ee <- eigen(Ka, symmetric = TRUE)
    ev <- ee$values
    cut <- max(length(ev), 1L) * .Machine$double.eps *
      max(1, max(abs(ev)))
    cut <- max(cut, noise)

    if(any(ev < -cut))
      .stop(paste0(
        "joint penalized information is indefinite; minimum eigenvalue = ",
        format(min(ev), scientific = TRUE),
        ". Use method = \"working\" only if that approximation is intended."
      ))

    null <- which(ev <= cut)
    nonestimable <- if(length(null)) {
      rowSums(ee$vectors[, null, drop = FALSE]^2) >
        sqrt(.Machine$double.eps)
    } else {
      rep(FALSE, length(active))
    }

    if(any(nonestimable)) {
      bad <- active[nonestimable]
      V[bad, ] <- V[, bad] <- NA_real_
    }

    keep <- active[!nonestimable]
    if(length(keep)) {
      Kg <- K[keep, keep, drop = FALSE]
      Kg <- 0.5 * (Kg + t(Kg))
      R <- tryCatch(chol(Kg), error = function(e) NULL)
      if(is.null(R))
        .stop("failed to factor the positive identifiable information block.")
      V[keep, keep] <- chol2inv(R)
    }

    dimnames(V) <- list(names, names)
    return(V)
  }

  .cor <- function(V)
  {
    d <- diag(V)
    s <- sqrt(d)
    C <- V / outer(s, s)
    ok <- is.finite(d) & d > 0
    diag(C)[ok] <- 1
    diag(C)[!ok] <- NA_real_
    dimnames(C) <- dimnames(V)
    return(C)
  }

  ## fitted state --------------------------------------------------------

  if(!inherits(object, "gamlss2"))
    .stop("'object' must inherit from class \"gamlss2\".")

  if(!is.logical(full) || length(full) != 1L || is.na(full))
    .stop("'full' must be TRUE or FALSE.")

  if(type == "coef")
    return(stats::coef(object, full = full, dropall = FALSE))

  if(inherits(object, "bamlss2"))
    .stop(paste0(
      "'vcov()' is a penalized-likelihood covariance calculation and is ",
      "not intended for posterior-simulation objects."
    ))

  family <- tryCatch(stats::family(object), error = function(e) NULL)
  if(is.null(family) || is.null(family$names) || is.null(family$links))
    .stop("cannot obtain the fitted gamlss2 family.")

  pn <- as.character(family$names)
  if(!length(pn))
    .stop("fitted family has no distributional parameter names.")

  if(is.null(object$fitted.values))
    .stop(paste0(
      "final additive predictors are not available. Models fitted with ",
      "'light = TRUE' must be refitted before joint covariance is computed."
    ))

  eta <- .eta_list(object$fitted.values, pn)
  n <- length(eta[[1L]])

  if(any(vapply(eta, length, integer(1L)) != n))
    .stop("fitted distributional predictors have incompatible lengths.")

  if(any(!vapply(eta, .finite, logical(1L))))
    .stop("fitted distributional predictors contain non-finite values.")

  y <- .response(n)
  prior <- .prior_weights(n)

  ## Parameter values at the fitted predictors. The family map is authoritative
  ## for custom and composite families; applying links one by one is not.
  if(!is.function(family$map2par))
    .stop("the fitted family has no predictor-to-parameter map.")
  par <- family$map2par(eta)
  if(!is.list(par))
    .stop("the fitted family returned an invalid parameter object.")
  for(j in seq_along(pn)) {
    z <- .get_named(par, pn[j], j)
    if(length(z) == 1L)
      z <- rep(z, n)
    if(length(z) != n || any(!is.finite(z)))
      .stop(paste0(
        "invalid fitted distributional parameter '", pn[j], "'."
      ))
    par[[pn[j]]] <- z
  }

  ## construct complete coefficient-linear predictor blocks
  blocks <- vector("list", length(pn))
  names(blocks) <- pn
  for(j in seq_along(pn))
    blocks[[j]] <- .parameter_block(pn[j], j, n)
  blocks <- blocks[!vapply(blocks, is.null, logical(1L))]

  if(!length(blocks))
    .stop("the fitted model contains no free coefficient-linear predictors.")

  bn <- names(blocks)
  pblock <- vapply(blocks, function(z) length(z$coefficients), integer(1L))
  starts <- cumsum(c(1L, head(pblock, -1L)))
  ends <- cumsum(pblock)
  index <- Map(seq.int, starts, ends)
  names(index) <- bn

  beta <- unlist(lapply(blocks, `[[`, "coefficients"),
    recursive = FALSE, use.names = FALSE)
  beta <- as.numeric(beta)

  coefficient.names <- unlist(lapply(blocks, `[[`, "names"),
    recursive = FALSE, use.names = FALSE)
  coefficient.names <- as.character(coefficient.names)

  linear <- unlist(Map(function(z, ii) {
    if(!length(z$linear))
      return(integer())
    ii[z$linear]
  }, blocks, index), use.names = FALSE)

  forced.na <- !is.finite(beta)
  fixed <- unlist(lapply(blocks, `[[`, "fixed"), use.names = FALSE)
  block.fixed <- vapply(blocks, function(z) all(z$fixed), logical(1L))
  beta0 <- beta
  beta0[forced.na] <- 0

  ## verify that the assembled coefficient order agrees with coef.gamlss2().
  ## This catches accidental parameter or smooth-block reordering early.
  has.fixed <- any(fixed)
  cf0 <- tryCatch(
    .coef_vector(stats::coef(object, full = TRUE)),
    error = function(e) NULL
  )

  if(!is.null(cf0)) {
    if(!has.fixed && length(cf0) != length(beta))
      .stop(paste0(
        "assembled full coefficient vector has length ", length(beta),
        " but coef(object, full = TRUE) has length ", length(cf0), "."
      ))

    if(length(cf0) == length(beta) &&
      !.same(as.numeric(cf0), beta))
      .stop(paste0(
        "assembled coefficient order does not agree with ",
        "coef(object, full = TRUE)."
      ))
  }

  ## use coef.gamlss2() only to adopt the public coefficient naming
  cf <- .canonical_coef(beta, TRUE, linear)
  if(length(cf) == length(beta) && .same(cf, beta)) {
    if(!is.null(names(cf)) && all(nzchar(names(cf))))
      coefficient.names <- names(cf)
  }

  names(beta) <- coefficient.names
  names(beta0) <- coefficient.names

  ## complete design and penalty matrices
  X <- lapply(blocks, `[[`, "X")
  P <- .block_diag(blocks, "P")
  dimnames(P) <- list(coefficient.names, coefficient.names)

  ## likelihood information ---------------------------------------------

  H <- matrix(0, length(beta), length(beta),
    dimnames = list(coefficient.names, coefficient.names))

  if(method == "joint") {
    if(is.null(family$score) || !is.list(family$score))
      .stop(paste0(
        "the fitted family has no completed score list. Linked-scale ",
        "scores created by complete_family() are required."
      ))

    for(a in seq_along(blocks)) {
      if(block.fixed[a])
        next
      na <- bn[a]
      ia <- index[[a]]
      Xa <- X[[a]]

      for(b in a:length(blocks)) {
        if(block.fixed[b])
          next
        nb <- bn[b]
        ib <- index[[b]]
        Xb <- X[[b]]

        h <- .observed_curvature(family, eta, y, na, nb, n)
        wh <- prior * h

        Hab <- crossprod(Xa, Xb * wh)

        H[ia, ib] <- Hab
        if(a != b)
          H[ib, ia] <- t(Hab)
      }
    }
  }

  if(method == "working") {
    for(a in seq_along(blocks)) {
      if(block.fixed[a])
        next
      na <- bn[a]
      ia <- index[[a]]
      Xa <- X[[a]]

      w <- .working_weights(family, par, y, eta[[na]], na, n)
      H[ia, ia] <- crossprod(Xa, Xa * (prior * w))
    }
  }

  if(method == "numeric") {
    if(!is.function(family$pdf))
      .stop("the fitted family has no density function for numeric Hessian.")

    base <- beta0
    active.numeric <- which(!forced.na & !fixed)

    ## eta(beta) is represented relative to the fitted predictor, which
    ## automatically retains all fitted offsets and fixed contributions
    nll <- function(b, ...)
    {
      theta <- base
      theta[active.numeric] <- b
      ee <- eta

      for(a in seq_along(blocks)) {
        nm <- bn[a]
        ii <- index[[a]]
        ee[[nm]] <- eta[[nm]] +
          drop(X[[a]] %*% (theta[ii] - base[ii]))
      }

      pp <- family$map2par(ee)
      ld <- family$pdf(par = pp, y = y, log = TRUE)
      ld <- as.numeric(ld)

      if(length(ld) == 1L && n == 1L)
        return(-prior * ld)

      if(length(ld) != n || any(!is.finite(ld)))
        return(Inf)

      return(-sum(prior * ld))
    }

    dots <- list(...)
    if(length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots)))))
      .stop("all control arguments in '...' must be named.")

    bad <- intersect(names(dots), c("par", "fn", "gr"))
    if(length(bad))
      .stop(paste0(
        "arguments in '...' must not override: ",
        paste(bad, collapse = ", "), "."
      ))

    Hn <- if(length(active.numeric)) {
      tryCatch(
        do.call(stats::optimHess,
          c(list(par = base[active.numeric], fn = nll), dots)),
        error = function(e) e
      )
    } else {
      matrix(0, 0L, 0L)
    }

    if(inherits(Hn, "error"))
      .stop(paste0(
        "numerical coefficient Hessian failed: ",
        conditionMessage(Hn)
      ))

    Hn <- .matrix(Hn)
    if(is.null(Hn) || any(dim(Hn) != length(active.numeric)) ||
      any(!is.finite(Hn)))
      .stop("numerical coefficient Hessian returned invalid values.")

    H <- matrix(0, length(beta), length(beta),
      dimnames = list(coefficient.names, coefficient.names))
    if(length(active.numeric))
      H[active.numeric, active.numeric] <- 0.5 * (Hn + t(Hn))
  }

  ## add complete fitted quadratic penalty
  K <- H + P
  K <- 0.5 * (K + t(K))

  ## The curvature calculations are numerical. Scale their round-off tolerance
  ## by the likelihood information, not by K: a very large smooth penalty must
  ## not erase a smaller but valid likelihood direction.
  noise <- 100 * .Machine$double.eps^(2 / 3) * max(1, max(abs(H)))
  V <- .invert(K, coefficient.names, forced.na, fixed, noise)

  ## full = FALSE is applied only after inversion, hence the ordinary
  ## linear coefficient covariance is marginalized over smooth coefficients
  keep <- if(full) seq_along(beta) else linear

  V <- V[keep, keep, drop = FALSE]
  b <- beta[keep]

  ## use public coefficient names where available
  bc <- .canonical_coef(beta, full, linear)
  if(length(bc) == length(b) && .same(bc, b) &&
    !is.null(names(bc)) && all(nzchar(names(bc)))) {
    names(b) <- names(bc)
    dimnames(V) <- list(names(bc), names(bc))
  } else {
    names(b) <- rownames(V)
  }

  ## return requested representation
  if(type == "se") {
    se <- sqrt(diag(V))
    names(se) <- rownames(V)
    return(se)
  }

  if(type == "cor")
    return(.cor(V))

  return(V)
}


## Shared covariance consumers ------------------------------------------

## Map fitted coefficient blocks to the public full coefficient ordering.
## The map is intentionally constructed from public coefficient names rather
## than from the old inference.R coefficient_structure() implementation.
vcov_coefficient_blocks <- function(object, coefficient.names = NULL)
{
  if(is.null(coefficient.names))
    coefficient.names <- names(coef(object, full = TRUE, dropall = FALSE))
  parameters <- object$family$names
  blocks <- indices <- setNames(vector("list", length(parameters)), parameters)

  for(j in parameters) {
    prefix <- paste0(j, ".p.")
    ind <- which(startsWith(coefficient.names, prefix))
    if(length(ind)) {
      local <- substring(coefficient.names[ind], nchar(prefix) + 1L)
      blocks[[j]][[1L]] <- list(type = "linear", label = NULL,
        names = local, index = ind, active = seq_along(ind), size = length(ind))
      indices[[j]] <- ind
    }

    for(term in names(object$fitted.specials[[j]])) {
      prefix <- paste0(j, ".s.", term, ".")
      ind <- which(startsWith(coefficient.names, prefix))
      if(!length(ind))
        next
      local <- substring(coefficient.names[ind], nchar(prefix) + 1L)
      blocks[[j]][[length(blocks[[j]]) + 1L]] <- list(
        type = "smooth", label = term, names = local, index = ind,
        active = seq_along(ind), size = length(ind))
      indices[[j]] <- c(indices[[j]], ind)
    }
  }
  list(blocks = blocks, indices = indices)
}

## A compact signature used only to reject stale reusable prediction caches.
vcov_state_signature <- function(object)
{
  list(
    coefficients = unclass(coef(object, full = TRUE, dropall = FALSE)),
    fitted.values = object$fitted.values,
    weights = object$weights,
    linear.penalties = lapply(object$fitted.linear, function(x) x$penalty),
    special.penalties = lapply(object$fitted.specials, function(x)
      lapply(x, function(z) list(lambdas = z$lambdas, penalty = z$penalty)))
  )
}

## A small common representation for consumers that require full joint
## covariance.  It is derived exclusively from the new vcov.gamlss2().
vcov_information <- function(object,
  method = c("joint", "working", "numeric"), ...)
{
  method <- match.arg(method)
  beta <- unclass(coef(object, full = TRUE, dropall = FALSE))
  V <- vcov.gamlss2(object, full = TRUE, method = method, ...)
  if(!identical(names(beta), rownames(V)))
    stop("internal covariance and coefficient order do not agree", call. = FALSE)
  structure(c(list(coefficients = beta, covariance = V,
    dimension = length(beta), method = method,
    state = vcov_state_signature(object)),
    vcov_coefficient_blocks(object, names(beta))),
    class = "gamlss2.vcov.information")
}

## Row-wise diag(A V A') without allowing unused NA covariance rows to
## contaminate otherwise estimable predictions.
vcov_variance <- function(A, information)
{
  V <- if(inherits(information, "gamlss2.vcov.information"))
    information$covariance else information
  A <- as.matrix(A)
  if(ncol(A) != nrow(V))
    stop("prediction design and covariance dimensions do not agree", call. = FALSE)
  answer <- numeric(nrow(A))
  if(!nrow(A) || !ncol(A))
    return(answer)

  bad <- !is.finite(diag(V))
  affected <- if(any(bad)) rowSums(abs(A[, bad, drop = FALSE])) > 0 else
    rep(FALSE, nrow(A))
  keep <- !bad & colSums(abs(A)) > 0
  if(any(keep)) {
    Vk <- V[keep, keep, drop = FALSE]
    if(any(!is.finite(Vk))) {
      answer[] <- NA_real_
      return(answer)
    }
    Ak <- A[, keep, drop = FALSE]
    answer <- rowSums((Ak %*% Vk) * Ak)
    scale <- pmax(1, rowSums(abs(Ak))^2 * max(abs(Vk)))
    tiny <- answer < 0 & answer > -sqrt(.Machine$double.eps) * scale
    answer[tiny] <- 0
    answer[answer < 0] <- NA_real_
  }
  answer[affected] <- NA_real_
  answer
}

## Gaussian draws from the full covariance. Fixed coefficient directions have
## exactly zero variance and remain at their fitted value.
vcov_draws <- function(information, R, antithetic = FALSE, center = FALSE)
{
  R <- as.integer(R)[1L]
  if(is.na(R) || R < 1L)
    stop("'R' must be a positive integer", call. = FALSE)
  V <- information$covariance
  beta <- information$coefficients
  if(any(!is.finite(diag(V))))
    stop("Gaussian coefficient draws are unavailable with non-estimable coefficients",
      call. = FALSE)
  if(any(diag(V) < 0))
    stop("Gaussian coefficient draws require nonnegative variances",
      call. = FALSE)
  active <- which(diag(V) > 0)
  D <- matrix(0, length(beta), R)
  if(length(active)) {
    C <- tryCatch(chol(V[active, active, drop = FALSE]),
      error = function(e) NULL)
    if(is.null(C))
      stop("Gaussian coefficient draws require positive definite covariance",
        call. = FALSE)
    if(isTRUE(antithetic)) {
      R0 <- ceiling(R / 2)
      Z <- matrix(rnorm(length(active) * R0), length(active), R0)
      Z <- cbind(Z, -Z)[, seq_len(R), drop = FALSE]
    } else {
      Z <- matrix(rnorm(length(active) * R), length(active), R)
    }
    D[active, ] <- t(C) %*% Z
    if(isTRUE(center))
      D[active, ] <- sweep(D[active, , drop = FALSE], 1L,
        rowMeans(D[active, , drop = FALSE]), "-")
  }
  sweep(D, 1L, beta, "+")
}

## Block-diagonal covariance stored by the individual fitting steps. This is
## the fast default for summaries and effect plots and deliberately contains
## no cross-term or cross-parameter covariance.
vcov_local <- function(object,
  type = c("vcov", "cor", "se", "coef"), full = FALSE)
{
  type <- match.arg(type)
  beta <- unclass(coef(object, full = TRUE, dropall = FALSE))
  if(type == "coef")
    return(coef(object, full = full, dropall = FALSE))
  map <- vcov_coefficient_blocks(object, names(beta))$blocks
  V <- matrix(0, length(beta), length(beta),
    dimnames = list(names(beta), names(beta)))

  for(j in object$family$names) {
    fixed <- object$control$fixed[[j]]
    fixed <- length(fixed) == 1L && !is.na(fixed) && isTRUE(as.logical(fixed))
    for(block in map[[j]]) {
      if(fixed || !length(block$index))
        next
      source <- if(is.null(block$label)) object$fitted.linear[[j]] else
        object$fitted.specials[[j]][[block$label]]
      B <- source$vcov
      if(is.null(B) || !identical(dim(as.matrix(B)), c(block$size, block$size))) {
        V[block$index, ] <- V[, block$index] <- NA_real_
      } else {
        V[block$index, block$index] <- as.matrix(B)
      }
    }
  }

  linear <- unlist(lapply(map, function(z) {
    k <- which(vapply(z, function(x) identical(x$type, "linear"), logical(1L)))
    if(!length(k)) integer() else z[[k[1L]]]$index
  }), use.names = FALSE)
  keep <- if(isTRUE(full)) seq_along(beta) else linear
  V <- V[keep, keep, drop = FALSE]
  if(type == "se") {
    answer <- sqrt(pmax(diag(V), 0))
    answer[is.na(diag(V)) | diag(V) < 0] <- NA_real_
    names(answer) <- rownames(V)
    return(answer)
  }
  if(type == "cor") {
    s <- sqrt(diag(V))
    V <- V / outer(s, s)
    diag(V)[is.finite(s) & s > 0] <- 1
  }
  V
}
