is_constant_col <- function(x, tol = 1e-8) {
  if(!is.numeric(x)) return(FALSE)
  x <- x[is.finite(x)]
  if(length(x) == 0L) return(TRUE)
  (max(x) - min(x)) < tol
}

## Normal scaling (center + per-column scaling), keeping (Intercept) unscaled.
elm_normal_scale <- function(X) {
  X <- as.matrix(X)

  p <- ncol(X)
  if(is.null(p) || p == 0L) stop("X must have at least one column")

  ## Columns to scale: those that are NOT constant (intercept columns stay untouched).
  j <- !apply(X, 2, is_constant_col)

  ## If nothing to scale, return identity.
  if(!any(j)) {
    return(function(X) {
      as.matrix(X)
    })
  }

  Xj <- X[, j, drop = FALSE]

  ## Center + scale (protect against NA/0/near-0 sd).
  mX <- colMeans(Xj, na.rm = TRUE)
  sdX <- apply(Xj, 2, sd, na.rm = TRUE)
  sdX[!is.finite(sdX) | sdX < 1e-12] <- 1

  function(X) {
    X <- as.matrix(X)

    if(ncol(X) != p) {
      stop(sprintf("expected %d columns, got %d", p, ncol(X)))
    }

    tmp <- X[, j, drop = FALSE]
    tmp <- sweep(tmp, 2, mX, "-")
    tmp <- sweep(tmp, 2, sdX, "/")

    X[, j] <- tmp
    X
  }
}


## Group scaling (center + QR-based orthonormalization), keeping (Intercept) unscaled.
elm_group_scale <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  cn_all <- colnames(X)

  has_int <- "(Intercept)" %in% cn_all
  cn <- cn_all[cn_all != "(Intercept)"]  # columns to scale

  ## If there is nothing to scale, return identity
  if(length(cn) == 0L) {
    function(X) {
      X <- as.matrix(X)
      X <- X[, cn_all, drop = FALSE]
      colnames(X) <- cn_all
      rownames(X) <- rownames(X)
      X
    }
  } else {
    mu <- setNames(colMeans(X[, cn, drop = FALSE]), cn)
    Xc <- sweep(X[, cn, drop = FALSE], 2, mu[cn], "-")

    decomp <- qr(Xc)

    R <- qr.R(decomp)
    Tmat <- qr.solve(R, diag(ncol(R))) * sqrt(n)

    function(X) {
      X <- as.matrix(X)

      Xs <- matrix(0, nrow = nrow(X), ncol = length(cn_all))
      colnames(Xs) <- cn_all
      rownames(Xs) <- rownames(X)

      if(has_int) {
        Xs[, "(Intercept)"] <- X[, "(Intercept)"]
      }

      Xn <- X[, cn, drop = FALSE]
      Xc <- sweep(Xn, 2, mu[cn], "-")
      Xn <- Xc %*% Tmat

      Xs[, cn] <- Xn
      Xs
    }
  }
}

elm_sample_weights <- function(Z, k, a = "tanh",
  target_sd = 1.0, slope_sd_1d = 1.0, q_range = c(0.05, 0.95),
  min_sd = 1e-12)
{
  Z <- as.matrix(Z)
  n_weights <- ncol(Z)
  if(n_weights < 2L) stop("Z must have at least intercept + 1 column")

  p_in <- n_weights - 1L
  fan_in <- max(1L, p_in)

  ## activation-aware base scale
  base_sd <- switch(a,
    "relu"      = sqrt(2 / fan_in),         # He
    "logistic"  = sqrt(2 / (fan_in + 1)),   # Xavier-like
    "tanh"      = sqrt(2 / (fan_in + 1)),   # Xavier-like
    "identity"  = 1 / sqrt(fan_in),
    sqrt(2 / (fan_in + 1))
  )

  Zin <- Z[, -1L, drop = FALSE]  # non-intercept inputs

  if(p_in == 1L) {
    ## 1D: place transitions across the observed x-range
    x1 <- drop(Zin)
    x1 <- x1[is.finite(x1)]
    if(length(x1) < 5L) stop("not enough finite data for 1D init")

    q <- stats::quantile(x1, probs = q_range, na.rm = TRUE, names = FALSE)
    if(!is.finite(q[1]) || !is.finite(q[2]) || abs(q[2] - q[1]) < min_sd) {
      q <- range(x1)
    }
    r <- q[2] - q[1]
    if(!is.finite(r) || r < min_sd) r <- 1

    ## slope scale adapted to data
    sx <- stats::sd(x1, na.rm = TRUE)
    if(!is.finite(sx) || sx < min_sd) sx <- 1
    slope_sd <- slope_sd_1d / sx

    thr <- stats::runif(k, min = q[1], max = q[2])

    W <- sapply(seq_len(k), function(jj) {
      w <- numeric(n_weights)

      a1 <- stats::rnorm(1, mean = 0, sd = slope_sd)
      if(!is.finite(a1) || abs(a1) < 1e-6) a1 <- sign(a1 + 1e-12) * 1e-6

      b <- -a1 * thr[jj]  # put transition at thr
      w[1L] <- b
      w[2L] <- a1
      w
    })

  } else {
    ## multivariate: sample -> rescale to target_sd -> bias via median centering
    W <- sapply(seq_len(k), function(jj) {
      w_in <- stats::rnorm(fan_in, mean = 0, sd = base_sd)

      u <- drop(Zin %*% w_in)
      su <- stats::sd(u, na.rm = TRUE)
      if(!is.finite(su) || su < min_sd) su <- 1
      w_in <- w_in * (target_sd / su)

      b <- -stats::median(drop(Zin %*% w_in), na.rm = TRUE)

      w <- numeric(n_weights)
      w[1L] <- b
      w[-1L] <- w_in
      w
    }, simplify = "matrix")
    dim(W) <- c(n_weights, k)
  }

  ## keep weights as matrix n_weights x k
  W
}

elm <- function(x, k = 50, a = "tanh", ...)
{
  call <- match.call()

  xn <- substitute(x)

  formula <- try(as.formula(xn), silent = TRUE)

  if(!inherits(formula, "try-error")) {
    xn <- all.vars(formula)
    environment(formula) <- environment(x)
  } else {
    xn <- all.vars(xn)
    formula <- NULL
  }

  ## List for setting up the special model term. 
  st <- list()
  st$control <- list(...)
  if(is.null(st$control$criterion))
    st$control$criterion <- "bic"
  if(is.null(st$control$scale))
    st$control$scale <- TRUE
  st$term <- xn 
  st$label <- gsub(" ", "", paste0("elm(", as.character(deparse(call[[2]])), ")"))

  if(!is.null(formula)) {
    st$Z <- model.matrix(formula, na.action = na.pass)
  } else {
    if(is.factor(x)) {
      st$is_factor <- TRUE
      st$lev <- levels(x)
      st$is_ordered <- inherits(x, "ordered")
      if(st$is_ordered) {
        x <- ordered(as.character(x), levels = st$lev)
      } else {
        x <- factor(as.character(x), levels = st$lev)
      }
      st$Z <- model.matrix(~x, na.action = na.pass)
      colnames(st$Z) <- gsub("x", "", colnames(st$Z))
    } else {
      st$Z <- cbind(1, x)
      colnames(st$Z) <- c("(Intercept)", xn)
    }
  }

  st$colnames <- colnames(st$Z)
  st$formula <- formula

  if(isTRUE(st$is_factor) && st$control$scale) {
    st$scale_fun <- elm_group_scale(st$Z)
  }

  if(st$control$scale && is.null(st$scale_fun)) {
    st$scale_fun <- elm_normal_scale(st$Z)
  }

  if(!is.null(st$scale_fun)) {
    st$Z <- st$scale_fun(st$Z)
  }

  st$activation <- switch(a,
    "logistic" = function(x) plogis(pmax(pmin(x, 35), -35)),
    "tanh"     = function(x) tanh(pmax(pmin(x, 35), -35)),
    "relu"     = function(x) pmax(x, 0),
    "identity" = function(x) x,
    stop("Unknown activation '", a, "'. Use one of: logistic, tanh, relu, identity.")
  )

  st$n_weights <- ncol(st$Z)

  st$weights <- elm_sample_weights(st$Z, k = k, a = a)
  st$X <- st$activation(st$Z %*% st$weights)

#  st$X <- apply(st$weights, 2, FUN = function(w) {
#    st$activation(st$Z %*% w)
#  })

  ## Center hidden-layer design matrix (columns)
  st$X_center <- colMeans(st$X, na.rm = TRUE)
  st$X <- sweep(st$X, 2, st$X_center, "-")
  st$S <- list(diag(ncol(st$X)))
  st$rank <- qr(st$X)$rank
  st$pred_class <- "elm.fitted"
  st$ncol <- ncol(st$X)
  st$keep <- c("formula", "term", "weights", "scale_fun", "control",
    "is_factor", "colnames", "activation", "X_center", "ncol", "lev")

  ## Assign the "special" class and the new class "n".
  class(st) <- c("special", "elm", "X %*% b", "mcmc")

  return(st) 
}

special_fit.elm <- function(x, z, w, control, transfer, ...)
{
  rval <- smooth.construct_wfit(x, z, w, control, transfer, ...)

  ## Arguments needed for prediction.
  rval[x$keep] <- x[x$keep]

  class(rval) <- "elm.fitted"

  return(rval)
}

special_predict.elm.fitted <- function(x, data, se.fit = FALSE, samples = NULL, ...) 
{
  ## Build design matrix.
  if(!is.null(x$formula)) {
    tt <- terms(x$formula)
    mf <- model.frame(tt, data = data, na.action = na.pass)
    Z  <- model.matrix(tt, mf)
  } else {
    if(isTRUE(x$is_factor)) {
      vals <- as.character(data[[x$term]])
      bad <- which(!is.na(vals) & !(vals %in% x$lev))
      if(length(bad)) {
        print(head(unique(vals[bad]), 20))
        stop("Found values not in training levels")
      }
      lev <- x$lev
      if(is.null(lev))
        lev <- levels(data[[x$term]])
      if(isTRUE(x$is_ordered)) {
        data[[x$term]] <- ordered(as.character(data[[x$term]]), levels = lev)
      } else {
        data[[x$term]] <- factor(as.character(data[[x$term]]), levels = lev)
      }
      Z <- model.matrix(~data[[x$term]], na.action = na.pass)
      colnames(Z) <- gsub("data[[x$term]]", "", colnames(Z), fixed = TRUE)
    } else {
      Z <- cbind(1, data[[x$term]])
    }
  }

  ## Scale.
  if(!is.null(x$scale_fun)) {
    Z <- x$scale_fun(Z)
  }

  X <- apply(x$weights, 2, FUN = function(w) {
    x$activation(Z %*% w)
  })

  if(!is.null(x$X_center)) {
    X <- sweep(X, 2, x$X_center, "-")
  }

  if(!is.null(samples)) {
    fit <- apply(samples, 1, function(beta) {
      X %*% beta
    })
  } else {
    fit <- drop(X %*% x$coefficients)

    if(se.fit) {
      se <- rowSums((X %*% x$vcov) * X)
      se <- 2 * sqrt(se)
      fit <- cbind("fit" = fit, "lower" = fit - se, "upper" = fit + se)
    }
  }

  return(fit)
}

