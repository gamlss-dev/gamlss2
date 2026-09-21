## State used to validate the lazy inference cache.
gamlss2_inference_state <- function(object)
{
  linear.penalties <- lapply(object$fitted.linear, function(x) x$penalty)
  penalties <- lapply(object$fitted.specials, function(x) {
    lapply(x, function(z) list(lambdas = z$lambdas, penalty = z$penalty))
  })
  list(
    coefficients = unclass(coef(object, full = TRUE, dropall = FALSE)),
    fitted.values = object$fitted.values,
    weights = object$weights,
    linear.penalties = linear.penalties,
    penalties = penalties
  )
}

## Construct coefficient blocks directly from the fitted model structure.
## Qualified coefficient names are output labels only and are never parsed.
gamlss2_coefficient_structure <- function(object)
{
  if(is.null(object$fitted.values))
    stop("joint covariance requires retained fitted values; refit with light = FALSE")

  y <- if(is.null(object$y)) {
    model.response(model.frame(object, keepresponse = TRUE))
  } else {
    object$y
  }
  X <- if(is.null(object[["x"]])) model.matrix(object) else object[["x"]]
  eta <- as.list(object$fitted.values)
  parameters <- object$family$names
  eta <- eta[parameters]
  n <- if(is.null(dim(y))) length(y) else nrow(y)
  if(!length(eta) || any(lengths(eta) != n))
    stop("invalid final predictor state in fitted object")

  blocks <- designs <- indices <- setNames(vector("list", length(parameters)), parameters)
  penalties <- list()
  map <- list()
  coefficients <- numeric()
  active.coefficients <- numeric()
  dimension <- global <- 0L

  is_fixed <- function(j) {
    fixed <- object$control$fixed
    if(is.null(fixed)) return(FALSE)
    if(is.null(names(fixed))) {
      k <- match(j, parameters)
      fixed <- if(length(fixed) >= k) fixed[[k]] else FALSE
    } else {
      fixed <- fixed[[j]]
    }
    length(fixed) == 1L && !is.na(fixed) && as.logical(fixed)
  }

  add_block <- function(j, type, term, B, beta, P = NULL) {
    if(!is.numeric(beta) || is.list(beta))
      stop("joint covariance: unsupported nonlinear coefficient block '",
        term, "' in parameter '", j, "'")
    local.names <- names(beta)
    beta <- as.numeric(beta)
    if(is.null(local.names)) local.names <- as.character(seq_along(beta))
    if(type == "smooth") {
      prefix <- paste0(j, ".s.", term, ".")
      while(any(k <- startsWith(local.names, prefix)))
        local.names[k] <- substring(local.names[k], nchar(prefix) + 1L)
    } else {
      prefix <- paste0(j, ".p.")
    }
    qualified <- paste0(prefix, local.names)
    names(beta) <- qualified

    if(!is.matrix(B) || nrow(B) != n || ncol(B) != length(beta))
      stop("invalid coefficient design for term '", term,
        "' in parameter '", j, "'")
    if(any(!is.finite(B)))
      stop("non-finite coefficient design for term '", term,
        "' in parameter '", j, "'")

    g <- global + seq_along(beta)
    global <<- global + length(beta)
    fixed <- is_fixed(j)
    keep <- !fixed & !is.na(beta)
    if(any(!is.finite(beta[keep])))
      stop("non-finite coefficient in parameter '", j, "'")
    ind <- integer(length(beta))
    if(any(keep)) {
      ind[keep] <- dimension + seq_len(sum(keep))
      dimension <<- dimension + sum(keep)
      designs[[j]] <<- cbind(designs[[j]], B[, keep, drop = FALSE])
      indices[[j]] <<- c(indices[[j]], ind[keep])
      active.coefficients <<- c(active.coefficients, beta[keep])
    }

    reason <- rep(NA_character_, length(beta))
    if(fixed) reason[] <- "fixed"
    reason[!fixed & is.na(beta)] <- "aliased"
    map[[length(map) + 1L]] <<- data.frame(
      parameter = rep(j, length(beta)),
      type = rep(type, length(beta)),
      term = rep(if(is.null(term)) NA_character_ else term, length(beta)),
      local = seq_along(beta), global = g, name = qualified,
      active = keep, design = ifelse(keep, ind, NA_integer_),
      penalty = rep(if(is.null(term)) NA_character_ else term, length(beta)),
      reason = reason, stringsAsFactors = FALSE
    )
    coefficients <<- c(coefficients, setNames(beta, qualified))
    blocks[[j]][[length(blocks[[j]]) + 1L]] <<- list(
      type = type, label = term, names = local.names, global = g,
      active = keep, index = ind[keep], size = length(beta), design = B
    )

    if(!is.null(P) && any(keep)) {
      P <- as.matrix(P)
      if(!identical(dim(P), c(length(beta), length(beta))) ||
          any(!is.finite(P)))
        stop("invalid penalty for term '", term,
          "' in parameter '", j, "'")
      penalties[[length(penalties) + 1L]] <<- list(
        index = ind[keep], P = P[keep, keep, drop = FALSE])
    }
  }

  for(j in parameters) {
    beta <- object$coefficients[[j]]
    if(length(beta)) {
      B <- X[, names(beta), drop = FALSE]
      P <- NULL
      if(isTRUE(object$control$ridge)) {
        penalized <- as.numeric(names(beta) != "(Intercept)")
        lambda <- object$fitted.linear[[j]]$penalty
        if(is.null(lambda) && any(penalized))
          stop("final linear ridge penalty is unavailable; refit the model")
        if(is.null(lambda)) lambda <- 0
        P <- diag(penalized * lambda, length(beta))
      }
      add_block(j, "linear", NULL, B, beta, P)
    }

    for(term in names(object$fitted.specials[[j]])) {
      fitted <- object$fitted.specials[[j]][[term]]
      beta <- fitted$coefficients
      if(is.null(beta)) next
      special <- object$specials[[term]]
      if(is.null(special) || is.null(special$X))
        stop("joint covariance: unsupported term '", term,
          "' in parameter '", j, "'")
      B <- special$X
      if("(Intercept)" %in% names(beta) &&
          !"(Intercept)" %in% colnames(B))
        B <- cbind("(Intercept)" = 1, B)
      if(!is.null(names(beta)) && !is.null(colnames(B)) &&
          all(names(beta) %in% colnames(B)))
        B <- B[, names(beta), drop = FALSE]
      if(!is.null(special$binning))
        B <- B[special$binning$match.index, , drop = FALSE]

      P <- fitted$penalty
      if(is.null(P)) {
        S <- special[["S", exact = TRUE]]
        if(length(S)) {
          if(!is.list(S) || is.matrix(S)) S <- list(S)
          lambda <- fitted$lambdas
          if(is.null(lambda) || !length(lambda))
            stop("fitted smoothing parameters are unavailable for term '",
              term, "'")
          lambda <- rep(lambda, length.out = length(S))
          P <- matrix(0, length(beta), length(beta))
          for(k in seq_along(S))
            P <- P + lambda[k] * as.matrix(S[[k]])
        }
      }
      add_block(j, "smooth", term, B, beta, P)
    }
  }

  map <- if(length(map)) do.call("rbind", map) else data.frame()
  reference <- unclass(coef(object, full = TRUE, dropall = FALSE))
  if(!identical(names(coefficients), names(reference)) ||
      length(coefficients) != length(reference))
    stop("internal coefficient map does not match fitted coefficient order")

  P <- matrix(0, dimension, dimension)
  for(penalty in penalties)
    P[penalty$index, penalty$index] <-
      P[penalty$index, penalty$index, drop = FALSE] + penalty$P
  dimnames(P) <- list(names(active.coefficients), names(active.coefficients))

  list(y = y, eta = eta, n = n, blocks = blocks, designs = designs,
    indices = indices, penalty = P, map = map, coefficients = coefficients,
    active.coefficients = active.coefficients, dimension = dimension)
}

## Observation-wise negative Hessian on the linked predictor scale.
gamlss2_curvature <- function(family, eta, par, y, a, b,
  method = c("joint", "working"))
{
  method <- match.arg(method)
  n <- length(eta[[1L]])
  if(method == "working") {
    if(a != b) return(numeric(n))
    h <- family$hessian[[a]]
    if(!is.function(h))
      stop("working Hessian is unavailable for parameter '", a, "'")
    return(as.numeric(h(par = par, y = y)))
  }

  score <- family$score[[a]]
  if(!is.function(score))
    stop("linked score is unavailable for parameter '", a, "'")
  step <- .Machine$double.eps^(1/3) * pmax(1, abs(eta[[b]]))
  upper <- lower <- eta
  upper[[b]] <- eta[[b]] + step
  lower[[b]] <- eta[[b]] - step
  su <- as.numeric(score(par = family$map2par(upper), y = y))
  sl <- as.numeric(score(par = family$map2par(lower), y = y))
  if(length(su) != n || length(sl) != n)
    stop("family score must return one value per observation")
  numerical <- -(su - sl) / (2 * step)

  ## Some legacy cross callbacks contain Fisher/working approximations rather
  ## than the raw observed derivative. Use them only when their convention
  ## agrees with the linked-score derivative.
  if(a != b) {
    h <- family$hessian[[paste0(a, ":", b)]]
    if(is.function(h)) {
      analytic <- as.numeric(h(par = par, y = y))
      scale <- pmax(1, abs(numerical), abs(analytic))
      if(length(analytic) == n && all(is.finite(analytic)) &&
          all(is.finite(numerical)) &&
          max(abs(analytic - numerical) / scale) < 1e-5)
        return(analytic)
    }
  }
  numerical
}

## Symmetric factorization with explicit handling of unidentified directions.
gamlss2_factor_information <- function(A)
{
  p <- ncol(A)
  if(!p)
    return(list(factor = A, covariance = A, rank = 0L,
      nonestimable = logical()))
  A <- 0.5 * (A + t(A))
  R <- tryCatch(chol(A), error = function(e) NULL)
  if(!is.null(R)) {
    tolerance <- max(p, 1L) * .Machine$double.eps
    if(rcond(R)^2 > tolerance)
      return(list(factor = R, covariance = NULL, rank = p,
        nonestimable = rep(FALSE, p)))
  }

  ev <- eigen(A, symmetric = TRUE)
  scale <- max(1, abs(ev$values))
  tolerance <- max(p, 1L) * .Machine$double.eps * scale
  if(any(ev$values < -tolerance)) {
    warning("joint penalized information is materially indefinite; covariance is unavailable",
      call. = FALSE)
    V <- matrix(NA_real_, p, p, dimnames = dimnames(A))
    return(list(factor = NULL, covariance = V,
      rank = sum(ev$values > tolerance), nonestimable = rep(TRUE, p)))
  }
  positive <- ev$values > tolerance
  if(!any(positive)) {
    warning("joint penalized information has rank zero", call. = FALSE)
    V <- matrix(NA_real_, p, p, dimnames = dimnames(A))
    return(list(factor = NULL, covariance = V,
      rank = 0L, nonestimable = rep(TRUE, p)))
  }

  V <- tcrossprod(sweep(ev$vectors[, positive, drop = FALSE], 2L,
    sqrt(ev$values[positive]), "/"))
  null <- ev$vectors[, !positive, drop = FALSE]
  nonestimable <- if(ncol(null))
    rowSums(null^2) > sqrt(.Machine$double.eps) else rep(FALSE, p)
  V[nonestimable, ] <- V[, nonestimable] <- NA_real_
  dimnames(V) <- dimnames(A)
  warning("joint penalized information is rank deficient; non-estimable coefficient variances are NA",
    call. = FALSE)
  list(factor = NULL, covariance = V, rank = sum(positive),
    nonestimable = nonestimable)
}

## Build and cache joint penalized information at the final predictors.
joint_information <- function(object, method = c("joint", "working", "numeric"),
  cache = TRUE, ...)
{
  method <- match.arg(method)
  control <- list(...)
  state <- gamlss2_inference_state(object)
  cache.env <- attr(object, ".inference.cache", exact = TRUE)
  if(isTRUE(cache) && is.environment(cache.env)) {
    cached <- cache.env[[method]]
    if(!is.null(cached) && identical(cached$control, control) &&
        identical(cached$state, state,
        num.eq = FALSE, single.NA = FALSE))
      return(cached$value)
  }

  z <- gamlss2_coefficient_structure(object)
  family <- object$family
  parameters <- family$names
  weights <- object$weights
  if(is.null(weights)) weights <- rep(1, z$n)
  if(length(weights) != z$n || any(!is.finite(weights)) || any(weights < 0))
    stop("joint covariance requires finite, nonnegative prior weights")

  if(method == "numeric") {
    beta <- z$active.coefficients
    if(!length(beta)) {
      information <- precision <- z$penalty
    } else {
      loglik <- function(theta) {
        eta <- z$eta
        for(j in parameters) {
          ind <- z$indices[[j]]
          if(length(ind))
            eta[[j]] <- eta[[j]] +
              drop(z$designs[[j]] %*% (theta[ind] - beta[ind]))
        }
        ll <- if(is.null(object$weights)) {
          family$log_likelihood(par = family$map2par(eta), y = z$y)
        } else {
          sum(family$pdf(par = family$map2par(eta), y = z$y,
            log = TRUE) * weights, na.rm = TRUE)
        }
        ll - 0.5 * drop(crossprod(theta, z$penalty %*% theta))
      }
      gradient <- function(theta) {
        eta <- z$eta
        for(j in parameters) {
          ind <- z$indices[[j]]
          if(length(ind))
            eta[[j]] <- eta[[j]] +
              drop(z$designs[[j]] %*% (theta[ind] - beta[ind]))
        }
        par <- family$map2par(eta)
        g <- numeric(length(theta))
        for(j in parameters) {
          ind <- z$indices[[j]]
          if(length(ind))
            g[ind] <- drop(crossprod(z$designs[[j]],
              weights * family$score[[j]](par = par, y = z$y)))
        }
        g - drop(z$penalty %*% theta)
      }
      H <- try(optimHess(beta, fn = loglik, gr = gradient,
        control = control), silent = TRUE)
      if(inherits(H, "try-error") || anyNA(H))
        H <- optimHess(beta, fn = loglik, control = control)
      information <- -0.5 * (H + t(H)) - z$penalty
      precision <- information + z$penalty
    }
  } else {
    information <- matrix(0, z$dimension, z$dimension)
    par <- family$map2par(z$eta)
    active <- parameters[lengths(z$indices) > 0L]
    for(ii in seq_along(active)) {
      a <- active[ii]
      ia <- z$indices[[a]]
      for(jj in seq.int(ii, length(active))) {
        b <- active[jj]
        ib <- z$indices[[b]]
        h <- gamlss2_curvature(family, z$eta, par, z$y, a, b, method)
        if(length(h) != z$n || any(!is.finite(h)))
          stop("non-finite family curvature for parameters '", a,
            "' and '", b, "'")
        block <- crossprod(z$designs[[a]], z$designs[[b]] * (weights * h))
        information[ia, ib] <- block
        if(a != b) information[ib, ia] <- t(block)
      }
    }
    information <- 0.5 * (information + t(information))
    precision <- information + z$penalty
  }
  precision <- 0.5 * (precision + t(precision))
  dimnames(information) <- dimnames(precision) <- dimnames(z$penalty)
  factor <- gamlss2_factor_information(precision)
  rval <- c(z[c("blocks", "designs", "indices", "map", "coefficients",
    "active.coefficients", "dimension")], list(
      information = information, penalty = z$penalty, precision = precision,
      factor = factor$factor, covariance = factor$covariance,
      rank = factor$rank, nonestimable = factor$nonestimable,
      method = method, state = state
    ))
  class(rval) <- "gamlss2.information"
  if(isTRUE(cache) && is.environment(cache.env))
    cache.env[[method]] <- list(state = state, control = control, value = rval)
  rval
}

gamlss2_information_vcov <- function(info)
{
  if(!is.null(info$covariance)) return(info$covariance)
  if(is.null(info$factor))
    return(matrix(NA_real_, info$dimension, info$dimension))
  V <- chol2inv(info$factor)
  dimnames(V) <- dimnames(info$precision)
  V
}

## Expand active covariance to the complete external coefficient ordering.
gamlss2_expand_vcov <- function(info)
{
  p <- nrow(info$map)
  V <- matrix(0, p, p, dimnames = list(info$map$name, info$map$name))
  active <- which(info$map$active)
  if(length(active)) V[active, active] <- gamlss2_information_vcov(info)
  aliased <- which(info$map$reason == "aliased")
  if(length(aliased)) V[aliased, ] <- V[, aliased] <- NA_real_
  V
}

## Prediction variances without constructing a dense inverse on the full-rank path.
gamlss2_information_variance <- function(A, info)
{
  n <- nrow(A)
  variance <- numeric(n)
  if(!n || !ncol(A)) return(variance)
  if(!is.null(info$factor)) {
    for(first in seq.int(1L, n, by = 1024L)) {
      rows <- seq.int(first, min(n, first + 1023L))
      B <- backsolve(info$factor, t(A[rows, , drop = FALSE]), transpose = TRUE)
      variance[rows] <- colSums(B * B)
    }
  } else {
    V <- gamlss2_information_vcov(info)
    variance <- rowSums((A %*% V) * A)
  }
  variance
}

## Draw from N(beta, A^-1) using the precision factor directly.
gamlss2_information_draws <- function(info, R)
{
  if(is.null(info$factor) || info$rank != info$dimension)
    stop("Gaussian coefficient draws require full-rank positive definite information")
  Z <- matrix(rnorm(info$dimension * R), info$dimension, R)
  D <- backsolve(info$factor, Z)
  sweep(D, 1L, info$active.coefficients, "+")
}
