## Image model term using convolution, ReLU, spatial-pyramid pooling and a
## weighted ridge output layer.

## Normalize c(height, width[, channels]).
.im_image_dim <- function(image_dim, what = "dim")
{
  if(is.null(image_dim))
    return(NULL)
  if(!is.numeric(image_dim) || !length(image_dim) %in% c(2L, 3L) ||
      any(!is.finite(image_dim)) || any(image_dim < 1) ||
      any(image_dim != as.integer(image_dim))) {
    stop("im(): ", what,
      " must contain two or three positive integers: height, width[, channels]")
  }
  image_dim <- as.integer(image_dim)
  if(length(image_dim) == 2L)
    image_dim <- c(image_dim, 1L)
  image_dim
}

## Convert supported image representations to n x H x W x C.
## Supported inputs are a flattened n x p matrix, an n x H x W[ x C]
## array, or a list whose elements are equally-sized H x W[ x C] images.
.im_as_4d <- function(x, image_dim = NULL)
{
  if(is.null(image_dim))
    image_dim <- attr(x, "image_dim", exact = TRUE)
  image_dim <- .im_image_dim(image_dim)

  if(is.list(x) && !is.data.frame(x) && is.null(dim(x))) {
    if(!length(x))
      stop("im(): the image list is empty")
    first.dim <- dim(x[[1L]])
    if(is.null(first.dim) || !length(first.dim) %in% c(2L, 3L))
      stop("im(): every image in a list must be an H x W matrix or H x W x C array")
    first.dim <- as.integer(first.dim)
    if(length(first.dim) == 2L)
      first.dim <- c(first.dim, 1L)
    if(is.null(image_dim))
      image_dim <- first.dim
    if(!identical(first.dim, image_dim))
      stop("im(): dim does not match the dimensions of the images")

    out <- array(NA_real_, dim = c(length(x), image_dim))
    for(i in seq_along(x)) {
      di <- dim(x[[i]])
      if(is.null(di) || !length(di) %in% c(2L, 3L))
        stop("im(): every image in a list must be an H x W matrix or H x W x C array")
      di <- as.integer(di)
      if(length(di) == 2L)
        di <- c(di, 1L)
      if(!identical(di, image_dim))
        stop("im(): all images in a list must have the same dimensions")
      if(!is.numeric(x[[i]]) && !is.logical(x[[i]]))
        stop("im(): image values must be numeric")
      out[i, , , ] <- array(as.numeric(x[[i]]), dim = image_dim)
    }
    return(out)
  }

  dx <- dim(x)
  if(is.null(dx))
    stop("im(): image input must be a matrix, array, or list of images")
  if(!is.numeric(x) && !is.logical(x))
    stop("im(): image values must be numeric")

  ## n x p flattened images.
  if(length(dx) == 2L) {
    if(is.null(image_dim))
      stop("im(): for an n x p matrix, supply dim = c(height, width[, channels])")
    if(prod(image_dim) != ncol(x))
      stop("im(): prod(dim) must equal ncol(x)")

    ## Rows of x are observations. Transposing first makes all pixels of an
    ## observation contiguous before forming H x W x C x n.
    out <- array(as.numeric(t(x)), dim = c(image_dim, nrow(x)))
    return(aperm(out, c(4L, 1L, 2L, 3L)))
  }

  ## n x H x W and n x H x W x C arrays.
  if(length(dx) %in% c(3L, 4L)) {
    actual.dim <- as.integer(dx[-1L])
    if(length(actual.dim) == 2L)
      actual.dim <- c(actual.dim, 1L)
    if(!is.null(image_dim) && !identical(actual.dim, image_dim))
      stop("im(): dim does not match the image array dimensions")
    if(length(dx) == 3L)
      return(array(x, dim = c(dx, 1L)))
    return(x)
  }

  stop("im(): image input must have dimensions n x p, n x H x W, or n x H x W x C")
}

## Scale channels separately using training-data statistics.
.im_scale_fit <- function(X, scale = TRUE)
{
  C <- dim(X)[4L]
  if(!isTRUE(scale))
    return(list(center = rep(0, C), scale = rep(1, C)))

  center <- scale.value <- numeric(C)
  for(cc in seq_len(C)) {
    xx <- as.numeric(X[, , , cc])
    center[cc] <- mean(xx)
    scale.value[cc] <- stats::sd(xx)
  }
  scale.value[!is.finite(scale.value) | scale.value <= 0] <- 1
  list(center = center, scale = scale.value)
}

.im_scale_apply <- function(X, scaling)
{
  ## Early prototype compatibility.
  if(is.null(scaling$center)) {
    scaling$center <- scaling$mean
    scaling$scale <- scaling$sd
  }
  X <- sweep(X, 4L, scaling$center, "-")
  sweep(X, 4L, scaling$scale, "/")
}

## Same-size zero-padded convolution.
## X: n x H x W x C; K: kh x kw x C x F.
.im_conv_same <- function(X, K, bias, cache = FALSE)
{
  d <- dim(X)
  dk <- dim(K)
  n <- d[1L]
  H <- d[2L]
  W <- d[3L]
  C <- d[4L]
  kh <- dk[1L]
  kw <- dk[2L]
  F <- dk[4L]

  if(is.loaded("im_conv_same", PACKAGE = "gamlss2")) {
    rval <- .Call(C_im_conv_same, X, K, bias)
    return(if(isTRUE(cache)) rval else rval$value)
  }

  ph1 <- floor((kh - 1L) / 2L)
  ph2 <- kh - 1L - ph1
  pw1 <- floor((kw - 1L) / 2L)
  pw2 <- kw - 1L - pw1

  XP <- array(0, dim = c(n, H + ph1 + ph2, W + pw1 + pw2, C))
  XP[, (ph1 + 1L):(ph1 + H), (pw1 + 1L):(pw1 + W), ] <- X

  Z <- array(0, dim = c(n, H, W, F))
  for(f in seq_len(F)) {
    Z[, , , f] <- bias[f]
    for(cc in seq_len(C)) {
      for(a in seq_len(kh)) {
        rr <- a:(a + H - 1L)
        for(b in seq_len(kw)) {
          ss <- b:(b + W - 1L)
          Xslice <- array(XP[, rr, ss, cc, drop = FALSE], dim = c(n, H, W))
          Zf <- array(Z[, , , f, drop = FALSE], dim = c(n, H, W))
          Z[, , , f] <- Zf + K[a, b, cc, f] * Xslice
        }
      }
    }
  }

  if(isTRUE(cache)) list(value = Z, padded = XP) else Z
}

## Pixel regions for spatial-pyramid pooling.
.im_pool_regions <- function(H, W, pool)
{
  regions <- list()
  k <- 1L
  for(g in pool) {
    row.start <- floor((seq_len(g) - 1L) * H / g) + 1L
    row.end <- floor(seq_len(g) * H / g)
    col.start <- floor((seq_len(g) - 1L) * W / g) + 1L
    col.end <- floor(seq_len(g) * W / g)
    for(ir in seq_len(g)) {
      for(ic in seq_len(g)) {
        regions[[k]] <- list(
          rows = row.start[ir]:row.end[ir],
          cols = col.start[ic]:col.end[ic]
        )
        k <- k + 1L
      }
    }
  }
  regions
}

.im_pool <- function(A, pool)
{
  d <- dim(A)
  n <- d[1L]
  F <- d[4L]
  regions <- .im_pool_regions(d[2L], d[3L], pool)
  P <- matrix(0, nrow = n, ncol = length(regions) * F)
  k <- 1L
  for(region in regions) {
    for(f in seq_len(F)) {
      values <- A[, region$rows, region$cols, f, drop = FALSE]
      P[, k] <- rowMeans(matrix(values, nrow = n))
      k <- k + 1L
    }
  }
  P
}

## Backpropagate through spatial-pyramid mean pooling.
.im_pool_back <- function(dP, activation_dim, pool)
{
  n <- activation_dim[1L]
  H <- activation_dim[2L]
  W <- activation_dim[3L]
  F <- activation_dim[4L]
  regions <- .im_pool_regions(H, W, pool)
  dA <- array(0, dim = activation_dim)
  k <- 1L
  for(region in regions) {
    area <- length(region$rows) * length(region$cols)
    for(f in seq_len(F)) {
      add <- array(dP[, k] / area,
        dim = c(n, length(region$rows), length(region$cols)))
      current <- array(dA[, region$rows, region$cols, f, drop = FALSE],
        dim = c(n, length(region$rows), length(region$cols)))
      dA[, region$rows, region$cols, f] <- current + add
      k <- k + 1L
    }
  }
  dA
}

## Initialize convolutional weights without unexpectedly changing the caller's
## random-number stream when an explicit seed is supplied.
.im_init_params <- function(X, filters, kernel, seed = NULL)
{
  if(!is.null(seed)) {
    had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if(had.seed)
      old.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if(had.seed) {
        assign(".Random.seed", old.seed, envir = .GlobalEnv)
      } else if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    })
    set.seed(seed)
  }

  C <- dim(X)[4L]
  kh <- kernel[1L]
  kw <- kernel[2L]
  F <- filters
  ## He initialization for ReLU.
  sdk <- sqrt(2 / (kh * kw * C))
  list(
    K = array(stats::rnorm(kh * kw * C * F, sd = sdk),
      dim = c(kh, kw, C, F)),
    b = rep(0, F),
    v = numeric(0),
    a = 0
  )
}

.im_zero_like <- function(p)
{
  lapply(p, function(z) z * 0)
}

.im_adam_step <- function(par, grad, state, lr, beta1, beta2, eps)
{
  if(is.null(state))
    state <- list(m = .im_zero_like(par), v = .im_zero_like(par), t = 0L)

  state$t <- state$t + 1L
  tt <- state$t
  for(nm in names(par)) {
    state$m[[nm]] <- beta1 * state$m[[nm]] + (1 - beta1) * grad[[nm]]
    state$v[[nm]] <- beta2 * state$v[[nm]] + (1 - beta2) * grad[[nm]]^2
    mh <- state$m[[nm]] / (1 - beta1^tt)
    vh <- state$v[[nm]] / (1 - beta2^tt)
    par[[nm]] <- par[[nm]] - lr * mh / (sqrt(vh) + eps)
  }
  list(par = par, state = state)
}

## Exact weighted ridge fit of the output layer. The intercept is not
## penalized; centering P makes this both stable and inexpensive.
.im_output_fit <- function(P, z, w, decay, edf = FALSE)
{
  sw <- sum(w)
  pbar <- drop(crossprod(w, P) / sw)
  zbar <- sum(w * z) / sw
  Pc <- sweep(P, 2L, pbar, "-")
  zw <- z - zbar
  B <- crossprod(Pc, Pc * (w / sw))
  rhs <- drop(crossprod(Pc, w * zw / sw))

  A <- B
  diag(A) <- diag(A) + decay
  R <- tryCatch(chol(A), error = function(e) NULL)
  if(!is.null(R)) {
    v <- backsolve(R, forwardsolve(t(R), rhs))
  } else {
    ee <- eigen(B, symmetric = TRUE)
    tol <- max(1, max(abs(ee$values))) * 1e-10
    denom <- ee$values + decay
    keep <- denom > tol
    v <- numeric(ncol(P))
    if(any(keep)) {
      rotated <- drop(crossprod(ee$vectors, rhs))
      v <- drop(ee$vectors[, keep, drop = FALSE] %*%
        (rotated[keep] / denom[keep]))
    }
  }
  a <- zbar - sum(pbar * v)

  rval <- list(a = a, v = v)
  if(isTRUE(edf)) {
    ev <- pmax(eigen(B, symmetric = TRUE, only.values = TRUE)$values, 0)
    tol <- max(1, max(ev)) * 1e-10
    rval$edf <- if(decay > 0) sum(ev / (ev + decay)) else sum(ev > tol)
  }
  rval
}

.im_check_control <- function(image_dim, filters, kernel, pool, epochs,
  learning_rate, decay, beta1, beta2, eps, clip)
{
  scalar.integer <- function(value, name, lower = 1L) {
    if(length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
        value < lower || value != as.integer(value))
      stop("im(): ", name, " must be an integer >= ", lower)
    as.integer(value)
  }

  filters <- scalar.integer(filters, "filters")
  epochs <- scalar.integer(epochs, "epochs", lower = 0L)
  if(length(kernel) == 1L)
    kernel <- rep(kernel, 2L)
  if(length(kernel) != 2L)
    stop("im(): kernel must be one integer or c(height, width)")
  kernel <- vapply(kernel, scalar.integer, integer(1L), name = "kernel")

  if(!is.numeric(pool) || !length(pool) || any(!is.finite(pool)) ||
      any(pool < 1) || any(pool != as.integer(pool)))
    stop("im(): pool must contain positive integers")
  pool <- sort(unique(as.integer(pool)))
  if(any(pool > min(image_dim[1:2])))
    stop("im(): pooling grid levels cannot exceed image height or width")

  positive <- c(learning_rate = learning_rate, eps = eps)
  if(length(learning_rate) != 1L || length(eps) != 1L ||
      any(!is.finite(positive)) || any(positive <= 0))
    stop("im(): learning_rate and eps must be positive")
  if(length(decay) != 1L || !is.finite(decay) || decay < 0)
    stop("im(): decay must be a non-negative number")
  if(length(beta1) != 1L || length(beta2) != 1L ||
      !is.finite(beta1) || !is.finite(beta2) ||
      beta1 < 0 || beta1 >= 1 || beta2 < 0 || beta2 >= 1)
    stop("im(): beta1 and beta2 must be in [0, 1)")
  if(!is.null(clip) && (length(clip) != 1L || !is.finite(clip) || clip <= 0))
    stop("im(): clip must be NULL or a positive number")

  list(filters = filters, kernel = kernel, pool = pool, epochs = epochs)
}

.im_check_misc_control <- function(scale, seed, trace)
{
  if(length(scale) != 1L || !is.logical(scale) || is.na(scale))
    stop("im(): scale must be TRUE or FALSE")
  if(length(trace) != 1L || !is.logical(trace) || is.na(trace))
    stop("im(): trace must be TRUE or FALSE")
  if(!is.null(seed) && (length(seed) != 1L || !is.numeric(seed) ||
      !is.finite(seed)))
    stop("im(): seed must be NULL or one finite number")
  invisible(NULL)
}

## CNNfit(): weighted least-squares image update used by im().
CNNfit <- function(x, z, w = rep(1, length(z)), model = NULL,
  image_dim = NULL, filters = 4L, kernel = c(3L, 3L),
  pool = c(1L, 2L), epochs = 8L, learning_rate = 3e-3,
  decay = 1e-3, beta1 = 0.9, beta2 = 0.999, eps = 1e-8,
  clip = 5, scale = TRUE, seed = NULL, trace = FALSE)
{
  X <- .im_as_4d(x, image_dim)
  storage.mode(X) <- "double"
  z <- as.numeric(z)
  w <- pmax(as.numeric(w), 0)

  if(dim(X)[1L] != length(z))
    stop("CNNfit(): number of images must equal length(z)")
  if(length(w) != length(z))
    stop("CNNfit(): length(w) must equal length(z)")
  if(any(!is.finite(z)) || any(!is.finite(w)))
    stop("CNNfit(): z and w must be finite")
  if(sum(w) <= 0)
    stop("CNNfit(): positive total weight required")
  if(any(!is.finite(X)))
    stop("CNNfit(): image values must be finite")

  image_dim <- dim(X)[2:4]
  .im_check_misc_control(scale, seed, trace)
  checked <- .im_check_control(image_dim, filters, kernel, pool, epochs,
    learning_rate, decay, beta1, beta2, eps, clip)
  filters <- checked$filters
  kernel <- checked$kernel
  pool <- checked$pool
  epochs <- checked$epochs

  if(is.null(model)) {
    scaling <- .im_scale_fit(X, scale)
    XS <- .im_scale_apply(X, scaling)
    par <- .im_init_params(XS, filters, kernel, seed)
  } else {
    if(!inherits(model, "base_cnn"))
      stop("CNNfit(): model must be a fitted base_cnn object")
    ## Early prototype compatibility.
    scaling <- model$scaling
    if(is.null(scaling))
      scaling <- model$scale
    XS <- .im_scale_apply(X, scaling)
    par <- model$par

    if(!identical(as.integer(model$image_dim), as.integer(image_dim)))
      stop("CNNfit(): image dimensions changed between updates")
    if(!identical(as.integer(dim(par$K)[1:2]), kernel) ||
        dim(par$K)[3L] != image_dim[3L] || dim(par$K)[4L] != filters)
      stop("CNNfit(): convolution architecture changed between updates")
    if(!identical(as.integer(model$pool), pool))
      stop("CNNfit(): pooling architecture changed between updates")
  }

  n <- length(z)
  H <- image_dim[1L]
  W <- image_dim[2L]
  sw <- sum(w)
  ## RS iterations have new local targets.
  adam <- NULL
  best <- NULL
  best.objective <- Inf

  ## Epoch zero evaluates the warm start (or random feature map). Retaining
  ## the best epoch protects an RS update from a poor gradient step.
  for(ep in 0:epochs) {
    conv <- .im_conv_same(XS, par$K, par$b, cache = ep < epochs)
    if(ep < epochs) {
      Z <- conv$value
      XP <- conv$padded
    } else {
      Z <- conv
    }
    A <- pmax(Z, 0)
    P <- .im_pool(A, pool)
    output <- .im_output_fit(P, z, w, decay)
    par$a <- output$a
    par$v <- output$v
    fit <- drop(par$a + P %*% par$v)
    loss <- 0.5 * sum(w * (z - fit)^2) / sw
    objective <- loss + 0.5 * decay * (sum(par$K^2) + sum(par$v^2))

    if(objective < best.objective) {
      best.objective <- objective
      best <- list(par = par, fit = fit, P = P, loss = loss)
    }
    if(ep == epochs)
      break

    ## Gradient of the profiled objective. Treating the exactly optimized
    ## output coefficients as fixed here is the envelope-theorem gradient.
    dy <- w * (fit - z) / sw
    dP <- tcrossprod(dy, par$v)
    dA <- .im_pool_back(dP, dim(A), pool)
    dZ <- dA * (Z > 0)

    if(is.loaded("im_conv_gradient", PACKAGE = "gamlss2")) {
      gradient <- .Call(C_im_conv_gradient, XP, dZ, par$K,
        as.numeric(decay))
      gK <- gradient$K
      gb <- gradient$b
    } else {
      dk <- dim(par$K)
      kh <- dk[1L]
      kw <- dk[2L]
      C <- dk[3L]
      F <- dk[4L]
      gK <- array(0, dim = dk)
      gb <- numeric(F)
      for(f in seq_len(F)) {
        dZf <- array(dZ[, , , f, drop = FALSE], dim = c(n, H, W))
        gb[f] <- sum(dZf)
        for(cc in seq_len(C)) {
          for(a in seq_len(kh)) {
            rr <- a:(a + H - 1L)
            for(b in seq_len(kw)) {
              ss <- b:(b + W - 1L)
              Xslice <- array(XP[, rr, ss, cc, drop = FALSE], dim = c(n, H, W))
              gK[a, b, cc, f] <- sum(dZf * Xslice) +
                decay * par$K[a, b, cc, f]
            }
          }
        }
      }
    }
    if(!is.null(clip)) {
      gK <- pmax(-clip, pmin(clip, gK))
      gb <- pmax(-clip, pmin(clip, gb))
    }

    convolution.par <- list(K = par$K, b = par$b)
    updated <- .im_adam_step(convolution.par, list(K = gK, b = gb), adam,
      learning_rate, beta1, beta2, eps)
    par$K <- updated$par$K
    par$b <- updated$par$b
    adam <- updated$state

    if(isTRUE(trace))
      cat(sprintf("CNNfit epoch %d: loss = %.8f, objective = %.8f\n",
        ep + 1L, loss, objective))
  }

  par <- best$par
  output <- .im_output_fit(best$P, z, w, decay, edf = TRUE)
  ## Conditional output-layer EDF plus the learned convolutional parameters.
  ## Cap at the number of informative observations minus one so a model term
  ## cannot create negative residual degrees of freedom by construction.
  conv.df <- length(par$K) + length(par$b)
  edf <- min(conv.df + output$edf, max(0, sum(w > 0) - 1L))

  model.out <- list(
    par = par,
    scaling = scaling,
    image_dim = image_dim,
    filters = dim(par$K)[4L],
    kernel = dim(par$K)[1:2],
    pool = pool,
    decay = decay
  )
  class(model.out) <- "base_cnn"

  list(
    fitted.values = best$fit,
    model = model.out,
    loss = best$loss,
    objective = best.objective,
    edf = edf
  )
}

predict.base_cnn <- function(object, newdata, image_dim = NULL, ...)
{
  if(is.null(image_dim))
    image_dim <- object$image_dim
  X <- .im_as_4d(newdata, image_dim)
  storage.mode(X) <- "double"
  if(!identical(as.integer(dim(X)[2:4]), as.integer(object$image_dim)))
    stop("predict.base_cnn(): new images have different dimensions from the training images")
  if(any(!is.finite(X)))
    stop("predict.base_cnn(): image values must be finite")
  scaling <- object$scaling
  if(is.null(scaling))
    scaling <- object$scale
  X <- .im_scale_apply(X, scaling)
  Z <- .im_conv_same(X, object$par$K, object$par$b)
  P <- .im_pool(pmax(Z, 0), object$pool)
  drop(object$par$a + P %*% object$par$v)
}

## Image special constructor.
im <- function(x, dim = NULL, ...)
{
  expr <- substitute(x)
  if(!is.symbol(expr))
    stop("im(): x must be a single image variable, e.g. im(image, dim = c(12, 12))")
  term <- as.character(expr)
  term.label <- deparse1(expr, backtick = TRUE)

  ctr <- list(...)
  if(length(ctr) && (is.null(names(ctr)) || any(!nzchar(names(ctr)))))
    stop("im(): all control arguments in ... must be named")
  allowed <- c("filters", "kernel", "pool", "epochs", "warm_epochs",
    "learning_rate", "decay", "beta1", "beta2", "eps", "clip",
    "scale", "seed", "trace")
  unknown <- setdiff(names(ctr), allowed)
  if(length(unknown))
    stop("im(): unused control argument", if(length(unknown) > 1L) "s" else "",
      ": ", paste(unknown, collapse = ", "))

  if(is.null(ctr$filters))
    ctr$filters <- 4L
  if(is.null(ctr$kernel))
    ctr$kernel <- c(3L, 3L)
  if(is.null(ctr$epochs))
    ctr$epochs <- 8L
  if(is.null(ctr$warm_epochs))
    ctr$warm_epochs <- 3L
  if(is.null(ctr$learning_rate))
    ctr$learning_rate <- 3e-3
  if(is.null(ctr$decay))
    ctr$decay <- 1e-3
  if(is.null(ctr$beta1))
    ctr$beta1 <- 0.9
  if(is.null(ctr$beta2))
    ctr$beta2 <- 0.999
  if(is.null(ctr$eps))
    ctr$eps <- 1e-8
  if(!"clip" %in% names(ctr))
    ctr$clip <- 5
  if(is.null(ctr$scale))
    ctr$scale <- TRUE
  if(is.null(ctr$trace))
    ctr$trace <- FALSE

  X <- .im_as_4d(x, dim)
  image_dim <- base::dim(X)[2:4]
  .im_check_misc_control(ctr$scale, ctr$seed, ctr$trace)
  if(is.null(ctr$pool)) {
    max.level <- min(image_dim[1:2])
    ctr$pool <- unique(pmin(c(1L, 2L, 4L), max.level))
  }
  checked <- .im_check_control(image_dim, ctr$filters, ctr$kernel,
    ctr$pool, ctr$epochs, ctr$learning_rate, ctr$decay, ctr$beta1,
    ctr$beta2, ctr$eps, ctr$clip)
  ctr$filters <- checked$filters
  ctr$kernel <- checked$kernel
  ctr$pool <- checked$pool
  ctr$epochs <- checked$epochs
  ctr$warm_epochs <- .im_check_control(image_dim, ctr$filters, ctr$kernel,
    ctr$pool, ctr$warm_epochs, ctr$learning_rate, ctr$decay, ctr$beta1,
    ctr$beta2, ctr$eps, ctr$clip)$epochs

  st <- list(
    X = X,
    term = term,
    expr = expr,
    image_dim = image_dim,
    control = ctr,
    label = paste0("im(", term.label, ")")
  )
  class(st) <- c("special", "im")
  st
}

## RS backfitting update. CNN parameters are warm-started from the previous
## iteration; the output layer is always refitted to the new working response.
special_fit.im <- function(x, z, w, control, ...)
{
  transfer <- list(...)$transfer
  old.model <- if(is.null(transfer)) NULL else transfer$model
  cc <- x$control
  epochs <- if(is.null(old.model)) cc$epochs else cc$warm_epochs

  X <- x$X
  ## special_terms() may have evaluated the constructor on unique rows when
  ## global gamlss2 binning is enabled. Expand those images before fitting.
  if(!is.null(x$binning))
    X <- X[x$binning$match.index, , , , drop = FALSE]

  fit <- CNNfit(
    x = X,
    z = z,
    w = w,
    model = old.model,
    image_dim = x$image_dim,
    filters = cc$filters,
    kernel = cc$kernel,
    pool = cc$pool,
    epochs = epochs,
    learning_rate = cc$learning_rate,
    decay = cc$decay,
    beta1 = cc$beta1,
    beta2 = cc$beta2,
    eps = cc$eps,
    clip = cc$clip,
    scale = cc$scale,
    seed = cc$seed,
    trace = cc$trace
  )

  rval <- list(
    model = fit$model,
    fitted.values = fit$fitted.values,
    transfer = list(model = fit$model),
    edf = fit$edf,
    term = x$term,
    image_dim = x$image_dim,
    loss = fit$loss,
    objective = fit$objective
  )

  ## Center the additive effect; the ordinary model intercept absorbs the
  ## constant. Store the same shift for out-of-sample prediction.
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift
  class(rval) <- "im.fitted"
  rval
}

.im_newdata <- function(data, term, image_dim)
{
  if(!is.list(data) || is.null(data[[term]]))
    stop("im() prediction: newdata must contain image variable '", term, "'")
  .im_as_4d(data[[term]], image_dim)
}

special_predict.im.fitted <- function(x, data, se.fit = FALSE, ...)
{
  X <- .im_newdata(data, x$term, x$image_dim)
  p <- predict(x$model, newdata = X, image_dim = x$image_dim) - x$shift
  if(isTRUE(se.fit))
    p <- data.frame(fit = p)
  p
}

