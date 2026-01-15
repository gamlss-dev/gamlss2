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
  min_sd = 1e-12, orthogonalize = NULL)
{
  Z <- as.matrix(Z)
  n_weights <- ncol(Z)
  if(n_weights < 2L) stop("Z must have at least intercept + 1 column")

  p_in <- n_weights - 1L
  fan_in <- max(1L, p_in)

  orthogonalize <- if(is.null(orthogonalize)) {
    (p_in > 1L && k >= fan_in)
  } else {
    isTRUE(orthogonalize)
  }

  ## map activations to sensible defaults
  target_sd0 <- switch(a,
    "logistic"   = 0.8,
    "tanh"       = 1.2,
    "atan"       = 1.2,
    "softsign"   = 1.2,
    "relu"       = 1.6,
    "leaky_relu" = 1.6,
    "elu"        = 1.6,
    "softplus"   = 1.6,
    "gaussian"   = 0.6,
    "laplace"    = 0.8,
    "sine"       = 1.0,
    "identity"   = 1.0,
    1.0
  )

  if(isTRUE(all.equal(target_sd, 1.0))) {
    target_sd <- target_sd0
  }

  ## activation-aware base scale
  base_sd <- switch(a,
    "relu"       = sqrt(2 / fan_in),
    "leaky_relu" = sqrt(2 / fan_in),
    "elu"        = sqrt(2 / fan_in),
    "softplus"   = sqrt(2 / fan_in),
    "logistic"   = sqrt(2 / (fan_in + 1)),
    "tanh"       = sqrt(2 / (fan_in + 1)),
    "atan"       = sqrt(2 / (fan_in + 1)),
    "softsign"   = sqrt(2 / (fan_in + 1)),
    "sine"       = 1 / sqrt(fan_in),
    "gaussian"   = 1 / sqrt(fan_in),
    "laplace"    = 1 / sqrt(fan_in),
    "identity"   = 1 / sqrt(fan_in),
    sqrt(2 / (fan_in + 1))
  )

  Zin <- Z[, -1L, drop = FALSE]

  bias_from_u <- function(u) {
    qq <- stats::runif(1, min = 0.2, max = 0.8)
    -stats::quantile(u, probs = qq, na.rm = TRUE,
                     names = FALSE, type = 7)
  }

  if(p_in == 1L) {
    ## 1D
    x1 <- drop(Zin)
    x1 <- x1[is.finite(x1)]
    if(length(x1) < 5L) stop("not enough finite data for 1D init")

    q <- stats::quantile(x1, probs = q_range,
                          na.rm = TRUE, names = FALSE)
    if(!is.finite(q[1]) || !is.finite(q[2]) ||
       abs(q[2] - q[1]) < min_sd) {
      q <- range(x1)
    }

    sx_sd <- stats::sd(x1, na.rm = TRUE)
    sx_mad <- stats::mad(x1, constant = 1, na.rm = TRUE)
    sx <- sx_sd
    if(is.finite(sx_mad) && sx_mad > min_sd)
      sx <- 0.5 * sx_sd + 0.5 * sx_mad
    if(!is.finite(sx) || sx < min_sd) sx <- 1

    slope_sd <- slope_sd_1d / sx
    thr <- stats::runif(k, min = q[1], max = q[2])

    W <- sapply(seq_len(k), function(jj) {
      w <- numeric(n_weights)

      if(a %in% c("gaussian", "laplace")) {
        a1 <- abs(stats::rnorm(1, 0, slope_sd))
        a1 <- max(a1, 0.05)
        if(!is.finite(a1) || a1 < 1e-6) a1 <- 1e-6
        b <- -a1 * thr[jj]
      } else {
        a1 <- stats::rnorm(1, 0, slope_sd)
        if(!is.finite(a1) || abs(a1) < 1e-6)
          a1 <- sign(a1 + 1e-12) * 1e-6
        b <- -a1 * thr[jj]
      }

      w[1L] <- b
      w[2L] <- a1
      w
    })

  } else {
    ## multivariate
    W <- matrix(NA_real_, nrow = n_weights, ncol = k)
    jj <- 1L

    while(jj <= k) {

      m <- min(fan_in, k - jj + 1L)

      if(orthogonalize) {
        V <- matrix(stats::rnorm(fan_in * m), fan_in, m)
        decomp <- qr(V)
        Q <- qr.Q(decomp, complete = FALSE)
        V <- Q[, seq_len(m), drop = FALSE]
      } else {
        V <- matrix(stats::rnorm(fan_in * m), fan_in, m)
        nv <- sqrt(colSums(V^2))
        nv[!is.finite(nv) | nv < min_sd] <- 1
        V <- sweep(V, 2, nv, "/")
      }

      V <- V * (base_sd * sqrt(fan_in))

      for(ii in seq_len(m)) {

        w_in <- V[, ii]

        u <- drop(Zin %*% w_in)
        su <- stats::sd(u, na.rm = TRUE)
        if(!is.finite(su) || su < min_sd) su <- 1
        w_in <- w_in * (target_sd / su)

        u2 <- drop(Zin %*% w_in)
        b <- bias_from_u(u2)

        w <- numeric(n_weights)
        w[1L] <- b
        w[-1L] <- w_in

        W[, jj] <- w
        jj <- jj + 1L
      }
    }
  }

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
  st$control$termselect <- isTRUE(st$control$elastic)
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
    "logistic"   = function(x) plogis(pmax(pmin(x, 35), -35)),
    "tanh"       = function(x) tanh(pmax(pmin(x, 35), -35)),
    "relu"       = function(x) pmax(x, 0),
    "leaky_relu" = function(x) ifelse(x > 0, x, 0.01 * x),
    "elu" = function(x) {
      x0 <- pmax(pmin(x, 35), -35)
      ifelse(x0 > 0, x0, exp(x0) - 1)
    },
    "softplus"   = function(x) log1p(exp(pmax(pmin(x, 35), -35))),
    "atan"       = function(x) atan(x),
    "softsign"   = function(x) x / (1 + abs(x)),
    "gaussian"   = function(x) exp(-x^2),
    "laplace"    = function(x) exp(-abs(x)),
    "sine"       = function(x) sin(x),
    "identity"   = function(x) x,
    stop("Unknown activation '", a, "'.")
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

  X <- x$activation(Z %*% x$weights)

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

