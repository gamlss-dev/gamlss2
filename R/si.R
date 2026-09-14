## Smooth image term.
##
## The image coefficient surface is represented by a tensor product
## of marginal B-spline bases:
##
##   beta(h, w) = sum_j sum_k theta_jk Bh_j(h) Bw_k(w)
##
## and the image contribution is
##
##   f_i = sum_h sum_w x_i(h, w) beta(h, w).
##
## Hence the model matrix contains the projections of every image
## onto the tensor-product basis.
##
## For multiple channels, one coefficient surface is fitted per
## channel, with common smoothing parameters.

si <- function(x, dim = NULL, k = c(8L, 8L),
  degree = c(3L, 3L), m = c(2L, 2L), sp = NULL, ...)
{
  expr <- substitute(x)

  if(!is.symbol(expr))
    stop("si(): x must be a single image variable, e.g. ",
      "si(image, dim = c(28, 28))")

  term <- as.character(expr)
  term.label <- deparse1(expr, backtick = TRUE)

  ## Use the same image conventions as im().
  X4 <- .im_as_4d(x, dim)
  image_dim <- base::dim(X4)[2:4]

  rep2 <- function(x, name, lower = 1L) {
    if(length(x) == 1L)
      x <- rep(x, 2L)

    if(length(x) != 2L ||
        !is.numeric(x) ||
        any(!is.finite(x)) ||
        any(x < lower) ||
        any(x != as.integer(x))) {
      stop("si(): ", name, " must contain one or two integers >= ",
        lower)
    }

    as.integer(x)
  }

  k <- rep2(k, "k")
  degree <- rep2(degree, "degree", lower = 0L)
  m <- rep2(m, "m")

  if(any(k > image_dim[1:2]))
    stop("si(): k cannot exceed the corresponding image dimension")

  if(any(k < degree + 1L))
    stop("si(): each k must be at least degree + 1")

  if(any(m >= k))
    stop("si(): each difference-penalty order m must be smaller than k")

  ## A supplied sp fixes the two marginal smoothing parameters.
  if(!is.null(sp)) {
    if(length(sp) == 1L)
      sp <- rep(sp, 2L)

    if(length(sp) != 2L ||
        !is.numeric(sp) ||
        any(!is.finite(sp)) ||
        any(sp < 0)) {
      stop("si(): sp must be NULL, one non-negative number, ",
        "or two non-negative numbers")
    }
  }

  control <- list(...)

  sx <- list(
    term = term,
    label = paste0("si(", term.label, ")"),
    by = "NA",

    ## This is a two-dimensional smooth over the image domain.
    dim = 2L,

    image_dim = image_dim,
    k = k,
    degree = degree,
    m = m,

    ## Important: x itself is matrix-valued, but si() handles the
    ## image summation internally.
    sumConv = FALSE
  )

  if(length(control))
    sx$control <- control

  if(!is.null(sp))
    sx$sp <- sp

  class(sx) <- "si.smooth.spec"

  sx
}


## Extract the matrix-valued image variable from the data supplied by
## smoothCon() or PredictMat().
.si_data_image <- function(data, term)
{
  if(!is.list(data) || is.null(data[[term]]))
    stop("si(): data must contain image variable '", term, "'")

  data[[term]]
}


## Project n images onto Bh(h) %*% Bw(w).
##
## This deliberately avoids explicitly constructing the large
##
##   (H * W) x (kh * kw)
##
## Kronecker basis matrix. Algebraically this is exactly the same as
##
##   X %*% kronecker(Bw, Bh),
##
## but requires much less memory for large images.
.si_project <- function(x, image_dim, Bh, Bw)
{
  ## model.frame() may drop the dimensions of an AsIs matrix column and
  ## provide its values as one flattened vector. Restore the n x p matrix
  ## before handing it to the common image conversion routine.
  if(is.null(dim(x)) && (is.numeric(x) || is.logical(x))) {
    npix <- prod(image_dim)
    if(length(x) %% npix != 0L)
      stop("si(): flattened image data has incompatible length")
    x <- matrix(x, nrow = length(x) / npix)
  }

  X4 <- .im_as_4d(x, image_dim)
  storage.mode(X4) <- "double"

  n <- base::dim(X4)[1L]
  H <- image_dim[1L]
  W <- image_dim[2L]
  C <- image_dim[3L]

  kh <- ncol(Bh)
  kw <- ncol(Bw)

  Z <- matrix(0, nrow = n, ncol = C * kh * kw)

  for(cc in seq_len(C)) {
    ## n x H x W for one channel.
    Xi <- array(
      X4[, , , cc, drop = FALSE],
      dim = c(n, H, W)
    )

    ## First contract the H dimension:
    ##
    ##   n x W x H  ->  (n W) x H
    ##
    ## then multiply by Bh.
    Xh <- matrix(
      aperm(Xi, c(1L, 3L, 2L)),
      nrow = n * W,
      ncol = H
    )

    Xh <- Xh %*% Bh

    ## Xh now represents n x W x kh.
    Xh <- array(Xh, dim = c(n, W, kh))

    ## Contract the W dimension.
    Xw <- matrix(
      aperm(Xh, c(1L, 3L, 2L)),
      nrow = n * kh,
      ncol = W
    )

    Xw <- Xw %*% Bw

    ## n x kh x kw; flatten with kh varying fastest.
    Xw <- array(Xw, dim = c(n, kh, kw))
    Xw <- matrix(Xw, nrow = n)

    jj <- (cc - 1L) * kh * kw + seq_len(kh * kw)
    Z[, jj] <- Xw
  }

  ## Helpful column names for debugging / coefficient extraction.
  nm <- as.vector(
    outer(
      seq_len(kh),
      seq_len(kw),
      function(h, w) paste0("h", h, ".w", w)
    )
  )

  if(C == 1L) {
    colnames(Z) <- nm
  } else {
    colnames(Z) <- unlist(
      lapply(seq_len(C), function(cc)
        paste0("c", cc, ".", nm)),
      use.names = FALSE
    )
  }

  Z
}


## Construct the tensor-product image smooth.
smooth.construct.si.smooth.spec <- function(object, data, knots)
{
  x <- .si_data_image(data, object$term)

  H <- object$image_dim[1L]
  W <- object$image_dim[2L]
  C <- object$image_dim[3L]

  kh <- object$k[1L]
  kw <- object$k[2L]

  dh <- object$degree[1L]
  dw <- object$degree[2L]

  mh <- object$m[1L]
  mw <- object$m[2L]

  ## Marginal B-spline bases on the fixed image grid.
  Bh <- splines::bs(
    seq_len(H),
    df = kh,
    degree = dh,
    intercept = TRUE
  )

  Bw <- splines::bs(
    seq_len(W),
    df = kw,
    degree = dw,
    intercept = TRUE
  )

  ## Store these: prediction uses exactly the training image basis.
  object$Bh <- Bh
  object$Bw <- Bw

  ## Tensor-product image design matrix.
  object$X <- .si_project(
    x,
    image_dim = object$image_dim,
    Bh = Bh,
    Bw = Bw
  )

  ## Marginal P-spline difference penalties.
  Dh <- diff(diag(kh), differences = mh)
  Dw <- diff(diag(kw), differences = mw)

  Ph <- crossprod(Dh)
  Pw <- crossprod(Dw)

  ## Coefficients are ordered:
  ##
  ##   h1.w1, h2.w1, ..., hkh.w1,
  ##   h1.w2, ...
  ##
  ## which matches kronecker(Bw, Bh).
  Sh <- kronecker(diag(kw), Ph)
  Sw <- kronecker(Pw, diag(kh))

  ## For multichannel images fit one spatial coefficient surface per
  ## channel. The two smoothing parameters are shared across channels.
  if(C > 1L) {
    IC <- diag(C)
    Sh <- kronecker(IC, Sh)
    Sw <- kronecker(IC, Sw)
  }

  storage.mode(Sh) <- "double"
  storage.mode(Sw) <- "double"

  object$S <- list(Sh, Sw)

  object$bs.dim <- ncol(object$X)

  ## Rank of each marginal penalty.
  object$rank <- vapply(
    object$S,
    function(S) qr(S, tol = 1e-10)$rank,
    integer(1L)
  )

  ## Null space of the combined tensor-product penalty.
  St <- Sh + Sw
  object$null.space.dim <-
    ncol(St) - qr(St, tol = 1e-10)$rank

  ## Let smoothCon() impose the usual sum-to-zero constraint.
  ## In particular, do NOT set C to a zero-row matrix as lin() does.
  object$side.constrain <- TRUE
  object$plot.me <- FALSE

  ## Matrix-valued x has already been handled by .si_project().
  object$sumConv <- FALSE

  class(object) <- "si.effect"

  object
}


## Prediction matrix.
##
## mgcv::PredictMat() will subsequently apply exactly the same
## identifiability transformation that smoothCon() applied during fitting.
Predict.matrix.si.effect <- function(object, data)
{
  x <- .si_data_image(data, object$term)

  .si_project(
    x,
    image_dim = object$image_dim,
    Bh = object$Bh,
    Bw = object$Bw
  )
}
