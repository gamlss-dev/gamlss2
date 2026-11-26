## Selected working-covariance blocks without a combined training design.
## The large-system path applies D'WD + P term by term and solves one linear
## right-hand side at a time. Fitting, EDFs, and joint inference are unchanged.

## Matrix is optional and need not be attached to the search path. Use its
## generics explicitly when a stored design or penalty is an S4 matrix.
vcov_local_diag <- function(x)
  if(inherits(x, "Matrix")) Matrix::diag(x) else diag(x)

vcov_local_colmeans <- function(x)
  if(inherits(x, "Matrix")) Matrix::colMeans(x) else colMeans(x)

vcov_local_colsums <- function(x)
  if(inherits(x, "Matrix")) Matrix::colSums(x) else colSums(x)

vcov_local_crossprod <- function(x, y = NULL)
  if(inherits(x, "Matrix") || inherits(y, "Matrix")) Matrix::crossprod(x, y) else
    crossprod(x, y)

vcov_local_finite <- function(x)
{
  if(inherits(x, "Matrix"))
    return(!"x" %in% methods::slotNames(x) ||
      all(is.finite(methods::slot(x, "x"))))
  size <- length(x)
  if(!size) return(TRUE)
  for(first in seq.int(1, size, by = 1000000)) {
    ii <- seq.int(first, min(size, first + 999999))
    if(any(!is.finite(x[ii]))) return(FALSE)
  }
  TRUE
}

vcov_local_zero <- function(q)
  if(requireNamespace("Matrix", quietly = TRUE)) Matrix::Diagonal(q, 0) else
    matrix(0, q, q)

vcov_local_penalty <- function(P, term)
{
  if(!vcov_local_finite(P))
    stop("non-finite penalty for term '", term, "'.", call. = FALSE)
  q <- nrow(P)
  scale <- 1
  if(inherits(P, "Matrix")) {
    scale <- max(scale, Matrix::norm(P, "M"))
  } else for(first in seq.int(1, length(P), by = 1000000)) {
    ii <- seq.int(first, min(length(P), first + 999999))
    scale <- max(scale, abs(P[ii]))
  }
  tolerance <- 1e-9 * scale
  diagonal <- inherits(P, "diagonalMatrix") ||
    (requireNamespace("Matrix", quietly = TRUE) && Matrix::isDiagonal(P))
  if(diagonal) {
    valid <- min(vcov_local_diag(P)) >= -tolerance
    if(q > 128L && requireNamespace("Matrix", quietly = TRUE))
      P <- Matrix::Diagonal(q, vcov_local_diag(P))
  } else {
    if(is.matrix(P) && q > 128L && requireNamespace("Matrix", quietly = TRUE)) {
      nonzero <- 0
      for(first in seq.int(1, length(P), by = 1000000)) {
        ii <- seq.int(first, min(length(P), first + 999999))
        nonzero <- nonzero + sum(P[ii] != 0)
      }
      if(nonzero < 0.15 * length(P)) P <- Matrix::Matrix(P, sparse = TRUE)
    }
    P <- 0.5 * (P + if(inherits(P, "Matrix")) Matrix::t(P) else t(P))
    if(inherits(P, "sparseMatrix") && q > 128L) {
      valid <- tryCatch({
        suppressWarnings(Matrix::Cholesky(Matrix::forceSymmetric(P) +
          Matrix::Diagonal(q, 2 * tolerance), LDL = FALSE, perm = TRUE))
        TRUE
      }, error = function(e) FALSE)
    } else {
      valid <- min(eigen(as.matrix(P), symmetric = TRUE,
        only.values = TRUE)$values) >= -tolerance
    }
  }
  if(!valid)
    stop("penalty for term '", term, "' is not positive semidefinite.",
      call. = FALSE)
  P
}

vcov_local_design_matrix <- function(X, n, index = NULL, center = NULL)
{
  q <- ncol(X)
  if(is.null(center)) center <- numeric(q)
  list(X = X, n = n, q = q, index = index, center = center, generator = NULL)
}

## Preserve binned rows and centering as operators, rather than expanding or
## copying the complete design. Sparse designs remain sparse.
vcov_local_design <- function(fit, source, b, n)
{
  q <- length(b)
  b0 <- as.numeric(b)
  b0[!is.finite(b0)] <- 0
  fv <- fit$fitted.values
  if(!is.null(fv)) fv <- as.numeric(fv)
  same <- function(a, z) length(a) == length(z) &&
    all(is.finite(a)) && all(is.finite(z)) &&
    max(abs(a - z), 0) <= 1e-6 * max(1, abs(a), abs(z))

  if(inherits(fit, "fk.fitted") && is.numeric(source$x) &&
      length(source$x) == n && length(fit$basis.center) == q) {
    generator <- function(ii) fk_basis(source$x[ii], fit$interior.knots,
      fit$boundary, fit$degree, fit$basis.center)$X.full
    ## Verify each bounded chunk, including affine coefficient dependence
    ## through the stored, coefficient-independent spline basis.
    for(first in seq.int(1L, n, by = max(1L, floor(1000000 / q)))) {
      ii <- seq.int(first, min(n, first + max(1L, floor(1000000 / q)) - 1L))
      X <- generator(ii)
      if(ncol(X) != q || !vcov_local_finite(X) ||
          (!is.null(fv) && !same(drop(X %*% b0), fv[ii])))
        return(NULL)
    }
    return(list(X = NULL, n = n, q = q, index = NULL,
      center = numeric(q), generator = generator))
  }

  candidates <- list()
  for(z in list(fit, source)) {
    if(!is.list(z)) next
    for(nm in c("X", "x", "model.matrix", "design"))
      if(!is.null(z[[nm]])) candidates[[length(candidates) + 1L]] <- z[[nm]]
    for(nm in c("smooth", "model", "fit")) {
      nested <- z[[nm]]
      if(!is.list(nested)) next
      for(field in c("X", "x", "model.matrix", "design"))
        if(!is.null(nested[[field]]))
          candidates[[length(candidates) + 1L]] <- nested[[field]]
    }
  }
  for(X in candidates) {
    if(!(is.matrix(X) && is.numeric(X)) && !inherits(X, "Matrix")) next
    if(ncol(X) != q || !vcov_local_finite(X)) next
    index <- NULL
    if(nrow(X) != n) {
      index <- source$binning$match.index
      if(length(index) != n || anyNA(index) ||
          any(index < 1L | index > nrow(X))) next
      index <- as.integer(index)
    }
    counts <- if(is.null(index)) NULL else tabulate(index, nbins = nrow(X))
    center <- if(is.null(counts)) as.numeric(vcov_local_colmeans(X)) else
      as.numeric(vcov_local_crossprod(X, counts)) / n
    prediction <- as.numeric(X %*% b0)
    if(!is.null(index)) prediction <- prediction[index]
    if(inherits(fit, "ms.fitted")) {
      if(!is.null(fv) && !same(prediction - sum(center * b0), fv)) next
    } else if(is.null(fv) || same(prediction, fv)) {
      center[] <- 0
    } else if(!same(prediction - sum(center * b0), fv)) {
      next
    }
    ## A sparse copy is useful for large random-effect and local-support
    ## bases. Count nonzeros in bounded pieces before allocating it.
    if(is.matrix(X) && q >= 32L && length(X) >= 100000L &&
        requireNamespace("Matrix", quietly = TRUE)) {
      nonzero <- 0
      for(first in seq.int(1, length(X), by = 1000000)) {
        ii <- seq.int(first, min(length(X), first + 999999))
        nonzero <- nonzero + sum(X[ii] != 0)
      }
      if(nonzero < 0.15 * length(X)) X <- Matrix::Matrix(X, sparse = TRUE)
    }
    answer <- vcov_local_design_matrix(X, n, index, center)
    if(inherits(fit, "ms.fitted")) attr(answer, "ms.center") <-
      if(is.null(counts)) as.numeric(vcov_local_colmeans(X)) else
        as.numeric(vcov_local_crossprod(X, counts)) / n
    return(answer)
  }
  NULL
}

## Compact custom specials may expose only coefficient-linear prediction.
## Verify that contract on bounded chunks and regenerate rows on demand;
## never allocate an n-by-q Jacobian just to discover an unsupported term.
vcov_local_predict_design <- function(fit, b, data, n, cache.limit = 1000000L)
{
  if(is.null(data) || NROW(data) != n) return(NULL)
  q <- length(b)
  b0 <- b
  b0[!is.finite(b0)] <- 0
  predict.coefficients <- function(cf, ii) {
    candidate <- fit
    candidate$coefficients <- cf
    value <- tryCatch(special_predict(candidate,
      data = data[ii, , drop = FALSE], se.fit = FALSE), error = function(e) NULL)
    if(is.data.frame(value)) {
      value <- if("fit" %in% names(value)) value$fit else
        if(ncol(value) == 1L) value[[1L]] else NULL
    }
    if(is.matrix(value) && ncol(value) == 1L) value <- value[, 1L]
    if(inherits(value, "Matrix")) {
      if(nrow(value) != length(ii) || ncol(value) != 1L) return(NULL)
      value <- as.numeric(value)
    }
    if(!is.numeric(value) || length(value) != length(ii) ||
        any(!is.finite(value))) return(NULL)
    as.numeric(value)
  }
  jacobian <- function(ii, verify = FALSE) {
    f0 <- predict.coefficients(b0, ii)
    if(is.null(f0)) return(NULL)
    X <- matrix(0, length(ii), q)
    for(j in seq_len(q)) {
      h <- .Machine$double.eps^(1 / 3) * (1 + abs(b0[j]))
      upper <- lower <- b0
      upper[j] <- upper[j] + h
      lower[j] <- lower[j] - h
      a <- predict.coefficients(upper, ii)
      z <- predict.coefficients(lower, ii)
      if(is.null(a) || is.null(z)) return(NULL)
      X[, j] <- (a - z) / (2 * h)
    }
    if(!vcov_local_finite(X)) return(NULL)
    if(verify) {
      for(direction in list(rep(1, q), rep(c(-1, 1), length.out = q),
          seq_len(q) / q)) {
        direction <- direction / max(1, sqrt(sum(direction^2)))
        h <- 0.05 * (1 + max(abs(b0)))
        actual <- predict.coefficients(b0 + h * direction, ii)
        expected <- f0 + h * as.numeric(X %*% direction)
        if(is.null(actual) || max(abs(actual - expected)) >
            1e-6 * max(1, abs(actual), abs(expected))) return(NULL)
      }
    }
    list(X = X, prediction = f0)
  }
  cached <- if(n * as.double(q) <= cache.limit) matrix(0, n, q) else NULL
  chunk <- max(1L, min(50000L, floor(1000000 / q)))
  total <- numeric(q)
  smallest <- Inf
  largest <- -Inf
  scale <- 1
  for(first in seq.int(1L, n, by = chunk)) {
    ii <- seq.int(first, min(n, first + chunk - 1L))
    result <- jacobian(ii, verify = TRUE)
    if(is.null(result)) return(NULL)
    raw <- as.numeric(result$X %*% b0)
    target <- if(is.null(fit$fitted.values)) result$prediction else fit$fitted.values[ii]
    if(length(target) != length(ii) || any(!is.finite(target))) return(NULL)
    difference <- raw - target
    smallest <- min(smallest, difference)
    largest <- max(largest, difference)
    scale <- max(scale, abs(raw), abs(target))
    total <- total + colSums(result$X)
    if(!is.null(cached)) cached[ii, ] <- result$X
  }
  center <- total / n
  if(max(abs(c(smallest, largest))) <= 1e-6 * scale) {
    center[] <- 0
  } else if(max(abs(c(smallest, largest) - sum(center * b0))) > 1e-6 * scale) {
    return(NULL)
  }
  if(!is.null(cached)) return(vcov_local_design_matrix(cached, n, center = center))
  last.rows <- last.X <- NULL
  generator <- function(ii) {
    if(identical(ii, last.rows)) return(last.X)
    result <- jacobian(ii)
    if(is.null(result))
      stop("custom special no longer provides a finite coefficient-linear design.",
        call. = FALSE)
    last.rows <<- ii
    last.X <<- result$X
    last.X
  }
  list(X = NULL, n = n, q = q, index = NULL, center = center, generator = generator)
}

vcov_local_rows <- function(design, ii)
{
  X <- if(is.null(design$generator)) {
    jj <- if(is.null(design$index)) ii else design$index[ii]
    design$X[jj, , drop = FALSE]
  } else design$generator(ii)
  if(any(design$center != 0))
    X <- sweep(as.matrix(X), 2L, design$center, "-")
  X
}

vcov_local_expand <- function(block, b)
{
  answer <- numeric(block$design$q)
  answer[block$active] <- if(is.null(block$basis)) b else
    as.numeric(block$basis %*% b)
  answer
}

vcov_local_expand_abs <- function(block, b)
{
  answer <- numeric(block$design$q)
  answer[block$active] <- if(is.null(block$basis)) abs(b) else
    as.numeric(abs(block$basis) %*% abs(b))
  answer
}

vcov_local_reduce <- function(block, b)
{
  b <- b[block$active]
  if(is.null(block$basis)) b else as.numeric(crossprod(block$basis, b))
}

vcov_local_block <- function(design, b, P, C, covariance = NULL)
{
  active <- which(is.finite(b))
  basis <- NULL
  if(length(active) && nrow(C)) {
    C <- C[, active, drop = FALSE]
    C <- C[rowSums(abs(C)) > 0, , drop = FALSE]
    ## Coordinate constraints, including fk's padded zero, need no dense basis.
    coordinate <- if(nrow(C)) rowSums(C != 0) == 1L else logical()
    if(any(coordinate)) {
      keep <- colSums(abs(C[coordinate, , drop = FALSE])) == 0
      active <- active[keep]
      C <- C[!coordinate, keep, drop = FALSE]
    }
    if(nrow(C) && length(active)) {
      norm <- sqrt(rowSums(C^2))
      C <- C[norm > 0, , drop = FALSE] / norm[norm > 0]
      decomposition <- qr(t(C), tol = 1e-10)
      basis <- qr.Q(decomposition, complete = TRUE)
      basis <- basis[, seq_len(ncol(basis)) > decomposition$rank, drop = FALSE]
    }
  }
  size <- if(is.null(basis)) length(active) else ncol(basis)
  reduced.penalty <- NULL
  reference <- attr(P, "local.reference", exact = TRUE)
  if(size && !is.null(reference) && max(abs(as.matrix(P))) > 1e6) {
    ## Rotate stiff penalties into their range and common null space. This
    ## avoids subtracting enormous, nearly equal penalty contributions when
    ## an unpenalized direction is represented in the original basis.
    Q <- if(is.null(basis)) diag(length(active)) else basis
    R <- crossprod(Q, reference[active, active, drop = FALSE] %*% Q)
    ee <- eigen(0.5 * (R + t(R)), symmetric = TRUE)
    cut <- max(1L, nrow(R)) * .Machine$double.eps * max(1, abs(ee$values))
    positive <- ee$values > cut
    range <- ee$vectors[, positive, drop = FALSE]
    null <- ee$vectors[, !positive, drop = FALSE]
    if(ncol(range)) {
      Pr <- crossprod(range, crossprod(Q,
        as.matrix(P)[active, active, drop = FALSE] %*% Q) %*% range)
      spectrum <- eigen(0.5 * (Pr + t(Pr)), symmetric = TRUE)
      roundoff <- max(1L, nrow(Pr)) * .Machine$double.eps *
        max(1, abs(spectrum$values))
      if(any(spectrum$values < -roundoff))
        stop("local penalty is indefinite on its declared range.", call. = FALSE)
      basis <- Q %*% cbind(null, range %*% spectrum$vectors)
      reduced.penalty <- c(rep(0, ncol(null)), pmax(spectrum$values, 0))
    } else {
      basis <- Q
      reduced.penalty <- numeric(size)
    }
  }
  list(design = design, active = active, basis = basis, size = size,
    penalty = P, reduced.penalty = reduced.penalty, covariance = covariance)
}

vcov_local_controls <- function(control, p)
{
  defaults <- list(dense.limit = 128L, chunk.size = 50000L,
    tolerance = 1e-11, maxit = max(200L, min(10000L, 2L * p)))
  if(is.null(control)) return(defaults)
  if(!is.list(control) || is.null(names(control)) ||
      any(!names(control) %in% names(defaults)))
    stop("invalid 'local.vcov' controls.", call. = FALSE)
  for(name in names(control)) defaults[[name]] <- control[[name]]
  for(name in c("dense.limit", "chunk.size", "maxit")) {
    value <- defaults[[name]]
    lower <- if(name == "dense.limit") 0 else 1
    if(!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < lower || value != floor(value) || value > .Machine$integer.max)
      stop("invalid local covariance control '", name, "'.", call. = FALSE)
  }
  value <- defaults$tolerance
  if(!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < .Machine$double.eps || value >= 1)
    stop("invalid local covariance control 'tolerance'.", call. = FALSE)
  defaults
}

vcov_working_linear <- function(object, parameter, n, w, linear.matrix,
  special.source, special.matrix, special.penalty, special.constraints)
{
  b <- object$coefficients[[parameter]]
  p0 <- length(b)
  X <- linear.matrix(parameter, b, n)
  P <- vcov_local_zero(p0)
  if(isTRUE(object$control$ridge)) {
    lambda <- object$fitted.linear[[parameter]]$penalty
    if(length(lambda) != 1L || !is.finite(lambda) || lambda < 0)
      stop("invalid final linear ridge penalty.", call. = FALSE)
    values <- lambda * (names(b) != "(Intercept)")
    P <- if(requireNamespace("Matrix", quietly = TRUE))
      Matrix::Diagonal(p0, values) else diag(values, p0)
  }
  blocks <- list(vcov_local_block(vcov_local_design_matrix(X, n), b, P,
    matrix(0, 0L, p0), object$fitted.linear[[parameter]]$vcov))
  for(term in names(object$fitted.specials[[parameter]])) {
    fit <- object$fitted.specials[[parameter]][[term]]
    if(!is.list(fit) || is.null(fit$coefficients))
      stop("special term '", term, "' has no coefficient-linear representation.",
        call. = FALSE)
    bs <- fit$coefficients
    if(is.list(bs)) bs <- unlist(bs, use.names = TRUE)
    if(is.matrix(bs) && (ncol(bs) == 1L || nrow(bs) == 1L)) bs <- as.numeric(bs)
    if(!is.numeric(bs) || !is.null(dim(bs)) || !length(bs))
      stop("invalid coefficient vector for special term '", term, "'.",
        call. = FALSE)
    source <- special.source(parameter, term)
    D <- special.matrix(parameter, term, fit, source, bs, n, compact = TRUE)
    P <- special.penalty(parameter, term, fit, source, length(bs))
    C <- special.constraints(parameter, term, fit, D, bs)
    blocks[[length(blocks) + 1L]] <- vcov_local_block(D, bs, P, C, fit$vcov)
  }
  sizes <- vapply(blocks, `[[`, integer(1L), "size")
  ends <- cumsum(sizes)
  starts <- ends - sizes
  for(k in seq_along(blocks))
    blocks[[k]]$index <- starts[k] + seq_len(sizes[k])
  p <- sum(sizes)
  control <- vcov_local_controls(object$control$local.vcov, p)
  answer <- matrix(NA_real_, p0, p0, dimnames = list(names(b), names(b)))
  linear <- blocks[[1L]]$active
  if(!length(linear)) return(answer)
  result <- vcov_local_solve(blocks, w, seq_along(linear), control)
  answer[linear, linear] <- result$covariance
  attr(answer, "local.solver") <- result$details
  answer
}

vcov_local_solve <- function(blocks, w, linear, control)
{
  p <- sum(vapply(blocks, `[[`, integer(1L), "size"))
  n <- length(w)
  if(!vcov_local_finite(w) || any(w < 0))
    stop("invalid final local working weights.", call. = FALSE)
  dense.width <- vapply(blocks, function(x) {
    if(inherits(x$design$X, c("sparseMatrix", "diagonalMatrix")) &&
        !any(x$design$center != 0) && is.null(x$basis)) 1L else x$design$q
  }, integer(1L))
  chunk <- max(1L, min(control$chunk.size,
    floor(1000000 / max(1L, max(dense.width)))))
  rows <- function(first) seq.int(first, min(n, first + chunk - 1L))
  reduced.rows <- function(block, ii) {
    X <- vcov_local_rows(block$design, ii)[, block$active, drop = FALSE]
    if(!is.null(block$basis)) X <- X %*% block$basis
    X
  }
  penalty <- function(block, x) {
    if(!is.null(block$reduced.penalty)) return(block$reduced.penalty * x)
    vcov_local_reduce(block, as.numeric(block$penalty %*%
      vcov_local_expand(block, x)))
  }

  ## Small systems use a bounded crossproduct and selected triangular solves.
  ## This also provides rank-revealing handling for small unidentified models.
  if(p <= control$dense.limit) {
    K <- matrix(0, p, p)
    for(first in seq.int(1L, n, by = chunk)) {
      ii <- rows(first)
      D <- do.call(cbind, lapply(blocks, function(block)
        as.matrix(reduced.rows(block, ii))))
      D <- as.matrix(D) * sqrt(w[ii])
      K <- K + crossprod(D)
    }
    for(block in blocks) {
      if(!block$size) next
      for(j in seq_len(block$size)) {
        e <- numeric(block$size)
        e[j] <- 1
        K[block$index, block$index[j]] <-
          K[block$index, block$index[j]] + penalty(block, e)
      }
    }
    K <- 0.5 * (K + t(K))
    scale <- sqrt(pmax(diag(K), .Machine$double.xmin))
    A <- K / outer(scale, scale)
    R <- tryCatch(chol(A), error = function(e) NULL)
    if(!is.null(R) && rcond(R) > sqrt(p * .Machine$double.eps)) {
      rhs <- matrix(0, p, length(linear))
      rhs[cbind(linear, seq_along(linear))] <- 1 / scale[linear]
      solved <- backsolve(R, forwardsolve(t(R), rhs)) / scale
      V <- solved[linear, , drop = FALSE]
    } else {
      ee <- eigen(A, symmetric = TRUE)
      cut <- p * .Machine$double.eps * max(1, abs(ee$values))
      if(any(ee$values < -cut))
        stop("local penalized information is indefinite.", call. = FALSE)
      keep <- ee$values > cut
      root <- sweep(ee$vectors[linear, keep, drop = FALSE], 2L,
        sqrt(ee$values[keep]), "/")
      V <- tcrossprod(root) / outer(scale[linear], scale[linear])
      bad <- if(any(!keep)) rowSums(ee$vectors[linear, !keep, drop = FALSE]^2) >
        sqrt(.Machine$double.eps) else rep(FALSE, length(linear))
      V[bad, ] <- V[, bad] <- NA_real_
    }
    return(list(covariance = V, details = list(method = "selected Cholesky",
      dimension = p, observations = n, iterations = 0L, residual = 0)))
  }

  streamed <- n > 4 * chunk ||
    any(vapply(blocks, function(x) !is.null(x$design$generator), logical(1L)))
  multiply <- function(x, absolute = FALSE) {
    result <- numeric(p)
    full <- lapply(blocks, function(block)
      if(absolute) vcov_local_expand_abs(block, x[block$index]) else
        vcov_local_expand(block, x[block$index]))
    reduce <- function(block, value) {
      if(!absolute) return(vcov_local_reduce(block, value))
      value <- value[block$active]
      if(is.null(block$basis)) value else
        as.numeric(crossprod(abs(block$basis), value))
    }
    if(streamed || absolute) {
      for(first in seq.int(1L, n, by = chunk)) {
        ii <- rows(first)
        eta <- numeric(length(ii))
        for(k in seq_along(blocks)) {
          X <- vcov_local_rows(blocks[[k]]$design, ii)
          if(absolute) X <- abs(X)
          eta <- eta + as.numeric(X %*% full[[k]])
        }
        eta <- eta * w[ii]
        for(k in seq_along(blocks)) {
          block <- blocks[[k]]
          X <- vcov_local_rows(block$design, ii)
          if(absolute) X <- abs(X)
          result[block$index] <- result[block$index] + reduce(block,
            as.numeric(vcov_local_crossprod(X, eta)))
        }
      }
    } else {
      eta <- numeric(n)
      for(k in seq_along(blocks)) {
        d <- blocks[[k]]$design
        value <- as.numeric(d$X %*% full[[k]])
        if(!is.null(d$index)) value <- value[d$index]
        eta <- eta + value - sum(d$center * full[[k]])
      }
      eta <- eta * w
      for(block in blocks) {
        d <- block$design
        value <- eta
        if(!is.null(d$index)) {
          sums <- rowsum(value, d$index, reorder = TRUE)
          value <- numeric(nrow(d$X))
          value[as.integer(rownames(sums))] <- sums[, 1L]
        }
        value <- as.numeric(vcov_local_crossprod(d$X, value)) - d$center * sum(eta)
        result[block$index] <- vcov_local_reduce(block, value)
      }
    }
    for(k in seq_along(blocks)) {
      block <- blocks[[k]]
      value <- if(absolute) {
        ## For a PSD penalty, |P_ij| <= sqrt(P_ii P_jj). This
        ## supplies a conservative error bound with only a q-vector.
        if(!is.null(block$reduced.penalty)) {
          block$reduced.penalty * abs(x[block$index])
        } else {
          diagonal <- sqrt(pmax(vcov_local_diag(block$penalty), 0))
          reduce(block, diagonal * sum(diagonal * full[[k]]))
        }
      } else penalty(block, x[block$index])
      result[block$index] <- result[block$index] + value
    }
    result
  }

  ## The fitter's conditional covariance is a preconditioner only. It does
  ## not determine the answer and need not use exactly the final weights.
  for(k in seq_along(blocks)) {
    block <- blocks[[k]]
    B <- block$covariance
    if(!is.null(B) && identical(dim(B), c(block$design$q, block$design$q)) &&
        vcov_local_finite(B)) {
      if(length(block$active) != nrow(B))
        B <- B[block$active, block$active, drop = FALSE]
      if(!is.null(block$basis)) B <- vcov_local_crossprod(block$basis, B %*% block$basis)
      if(block$size && (any(vcov_local_diag(B) <= 0) ||
          (block$size <= 128L && is.null(tryCatch(chol(as.matrix(B)),
            error = function(e) NULL))))) B <- NULL
    } else B <- NULL
    if(is.null(B)) {
      diagonal <- numeric(block$size)
      for(first in seq.int(1L, n, by = chunk)) {
        ii <- rows(first)
        X <- reduced.rows(block, ii)
        diagonal <- diagonal + as.numeric(vcov_local_colsums(X^2 * w[ii]))
      }
      if(is.null(block$basis)) {
        diagonal <- diagonal + vcov_local_diag(block$penalty)[block$active]
      } else for(j in seq_len(block$size)) {
        e <- numeric(block$size)
        e[j] <- 1
        diagonal[j] <- diagonal[j] + penalty(block, e)[j]
      }
      diagonal.floor <- max(.Machine$double.xmin,
        .Machine$double.eps * max(diagonal, 0))
      B <- 1 / pmax(diagonal, diagonal.floor)
    }
    blocks[[k]]$preconditioner <- B
  }
  precondition <- function(r) {
    z <- numeric(p)
    for(block in blocks) {
      B <- block$preconditioner
      z[block$index] <- if(is.null(dim(B))) B * r[block$index] else
        as.numeric(B %*% r[block$index])
    }
    z
  }
  V <- matrix(0, length(linear), length(linear))
  iterations <- numeric(length(linear))
  residual <- numeric(length(linear))
  backward <- unscaled <- numeric(length(linear))
  for(column in seq_along(linear)) {
    rhs <- numeric(p)
    rhs[linear[column]] <- 1
    result <- vcov_local_pcg(multiply, precondition, rhs,
      control$tolerance, control$maxit,
      absolute = function(x) multiply(x, absolute = TRUE),
      roundoff = 8 * .Machine$double.eps * max(n, p))
    if(!result$converged)
      stop("local covariance solve did not converge (relative residual ",
        format(result$residual, scientific = TRUE), "); no dense fallback was used.",
        call. = FALSE)
    V[, column] <- result$solution[linear]
    iterations[column] <- result$iterations
    residual[column] <- result$residual
    unscaled[column] <- result$unscaled.residual
    backward[column] <- result$backward.error
  }
  scale <- sqrt(pmax(diag(V), .Machine$double.xmin))
  if(any(!is.finite(V)) || any(diag(V) <= 0) ||
      max(abs(V - t(V)) / outer(scale, scale)) > 1e-6)
    stop("local covariance solves failed their symmetry or positivity check.",
      call. = FALSE)
  list(covariance = 0.5 * (V + t(V)), details = list(method = "matrix-free PCG",
    dimension = p, observations = n, streamed = streamed,
    iterations = iterations, residual = residual,
    unscaled.residual = unscaled, backward.error = backward))
}

vcov_local_pcg <- function(multiply, precondition, rhs, tolerance, maxit,
  absolute, roundoff)
{
  x <- numeric(length(rhs))
  r <- rhs
  z <- precondition(r)
  direction <- z
  rho <- sum(r * z)
  norm <- sqrt(sum(rhs * z))
  unscaled.norm <- sqrt(sum(rhs^2))
  relative <- Inf
  check <- function(x, iteration) {
    r <- rhs - multiply(x)
    z <- precondition(r)
    energy <- sum(r * z)
    relative <- if(is.finite(energy) && energy >= 0) sqrt(energy) / norm else Inf
    unscaled <- sqrt(sum(r^2)) / unscaled.norm
    backward <- unscaled
    if(relative <= tolerance && unscaled > tolerance) {
      denominator <- abs(rhs) + absolute(x)
      backward <- max(abs(r) / pmax(denominator, .Machine$double.xmin))
    }
    list(solution = x, converged = is.finite(relative) && relative <= tolerance &&
      (unscaled <= tolerance || backward <= max(tolerance, roundoff)),
      iterations = iteration, residual = relative,
      unscaled.residual = unscaled, backward.error = backward, r = r, z = z)
  }
  for(iteration in seq_len(maxit)) {
    q <- multiply(direction)
    curvature <- sum(direction * q)
    if(!is.finite(curvature) || curvature <= 0 || !is.finite(rho) || rho <= 0) break
    alpha <- rho / curvature
    x <- x + alpha * direction
    r <- r - alpha * q
    z <- precondition(r)
    next.rho <- sum(r * z)
    relative <- if(is.finite(next.rho) && next.rho >= 0)
      sqrt(next.rho) / norm else Inf
    ## Check the true residual before accepting; restart periodically to
    ## prevent accumulated recurrence error from masquerading as convergence.
    restart <- iteration %% 50L == 0L || relative <= tolerance
    if(restart) {
      checked <- check(x, iteration)
      if(checked$converged) return(checked)
      r <- checked$r
      z <- checked$z
      next.rho <- sum(r * z)
    }
    direction <- if(restart) z else z + (next.rho / rho) * direction
    rho <- next.rho
  }
  check(x, iteration)
}
