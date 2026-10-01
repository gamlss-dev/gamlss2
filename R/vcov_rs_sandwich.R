## Post-fit sandwich covariance for the RS local REML fixed point.
## Case multipliers act on independent weighted rows. The local REML
## equations are evaluated at fixed coefficients without another smooth fit.
vcov_rs_sandwich <- function(object, beta, V, P, blocks, index, X, eta, y,
  prior, family, R, details = FALSE)
{
  if(!identical(object$control$optimizer, RS) || !is.null(object$jr) ||
      (!is.null(object$control$CG) &&
        !identical(object$control$CG, FALSE)))
    stop("RS sandwich covariance requires an ordinary RS fit.",
      call. = FALSE)
  if(any(!is.finite(diag(V))))
    stop("RS sandwich covariance requires estimable coefficients.",
      call. = FALSE)
  if(isTRUE(object$control$termselect) || isTRUE(object$control$ridge))
    stop("RS sandwich covariance does not support selected penalties.",
      call. = FALSE)
  if(!isTRUE(object$converged))
    stop("RS sandwich covariance requires a converged fit.",
      call. = FALSE)

  state <- vcov_unconditional(object, beta, V, P, blocks, index, R,
    state = TRUE)
  active <- state$active
  all.components <- state$components
  at.boundary <- vapply(all.components, function(z) {
    lambda <- object$fitted.specials[[z$parameter]][[z$term]]$
      lambdas[z$penalty]
    lambda == 1e-7 || lambda == 1e7
  }, logical(1L))
  boundary.components <- all.components[at.boundary]
  components <- all.components[!at.boundary]
  state$u <- state$u[, !at.boundary, drop = FALSE]
  m <- length(components)
  cn <- vapply(components, `[[`, character(1L), "name")
  rho <- if(m) log(vapply(components, function(z)
    object$fitted.specials[[z$parameter]][[z$term]]$
      lambdas[z$penalty], numeric(1L))) else numeric()
  names(rho) <- cn
  boundary.rho <- if(length(boundary.components))
    log(vapply(boundary.components, function(z)
      object$fitted.specials[[z$parameter]][[z$term]]$
        lambdas[z$penalty], numeric(1L))) else numeric()
  names(boundary.rho) <- vapply(boundary.components, `[[`,
    character(1L), "name")
  lookup <- integer(length(beta))
  lookup[active] <- seq_along(active)
  n <- length(prior)
  if(any(!is.finite(prior)) || any(prior < 0))
    stop("RS sandwich covariance requires finite nonnegative weights.",
      call. = FALSE)

  smooths <- list()
  for(a in seq_along(blocks)) {
    name <- names(blocks)[a]
    for(term in names(object$fitted.specials[[name]])) {
      ids <- which(vapply(components, function(z)
        identical(z$parameter, name) && identical(z$term, term),
        logical(1L)))
      boundary.ids <- which(vapply(boundary.components, function(z)
        identical(z$parameter, name) && identical(z$term, term),
        logical(1L)))
      if(!length(ids) && !length(boundary.ids)) next
      source <- object$specials[[term]]
      if(is.null(source) && !is.null(object$specials[[name]]))
        source <- object$specials[[name]][[term]]
      if(is.null(source) || !inherits(source, "mgcv.smooth") ||
          !is.null(source$special.wfit))
        stop("RS sandwich covariance requires ordinary mgcv smooths.",
          call. = FALSE)
      criterion <- source$control$criterion
      if(is.null(criterion)) criterion <- source$control$method
      if(is.null(criterion)) criterion <- object$control$criterion
      if(is.null(criterion))
        criterion <- if(isTRUE(object$control$logLik)) "aicc" else "reml"
      if(!identical(tolower(criterion), "reml") ||
          isTRUE(source$localML) || !is.null(source$sp))
        stop("RS sandwich covariance requires local REML smoothing estimates.",
          call. = FALSE)
      local <- which(blocks[[a]]$term == term)
      ii <- index[[a]][local]
      S <- lapply(source$S, as.matrix)
      lambda <- as.numeric(object$fitted.specials[[name]][[term]]$lambdas)
      rank <- if(length(source$null.space.dim) == 1L)
        length(ii) - source$null.space.dim else NULL
      basis <- if(length(S) > 1L) smooth.construct_reml(S, rank) else NULL
      if(length(S) > 1L && is.null(basis))
        stop("RS sandwich covariance requires an identifiable penalty range.",
          call. = FALSE)
      smooths[[length(smooths) + 1L]] <- list(
        parameter = name, index = ii, X = X[[a]][, local, drop = FALSE],
        S = S, lambda = lambda, ids = ids,
        boundary = vapply(boundary.ids, function(i)
          boundary.components[[i]]$penalty, integer(1L)),
        boundary.ids = boundary.ids, basis = basis,
        null = if(length(S) == 1L) {
          if(length(source$null.space.dim) == 1L) {
            source$null.space.dim
          } else {
            values <- eigen(S[[1L]], symmetric = TRUE,
              only.values = TRUE)$values
            sum(abs(values) <= sqrt(.Machine$double.eps) *
              max(1, max(abs(values))))
          }
        } else NULL)
    }
  }
  if(length(smooths)) {
    par <- family$map2par(eta)
    for(name in unique(vapply(smooths, `[[`, character(1L),
        "parameter"))) {
      ew <- .update(par = par, y = y, eta = eta[[name]],
        family = family, which = name)
      s <- as.numeric(family$score[[name]](par = par, y = y))
      working <- as.numeric(ew$weights) *
        (as.numeric(ew$eta) - eta[[name]])
      if(length(s) != n || length(working) != n ||
          any(!is.finite(working)) ||
          max(abs(s - working)) > 1e-6 * max(1, abs(s)))
        stop("RS working update does not match the linked coefficient score.",
          call. = FALSE)
    }
  }

  ## The local response is the final predictor with this smooth removed,
  ## plus the linked-scale score divided by its working curvature.
  equations <- function(theta, r, cases = FALSE, case.weights = prior)
  {
    ee <- eta
    for(a in seq_along(blocks)) {
      ii <- index[[a]]
      ee[[names(blocks)[a]]] <- eta[[names(blocks)[a]]] +
        drop(X[[a]] %*% (theta[ii] - beta[ii]))
    }
    par <- family$map2par(ee)
    out <- numeric(m)
    C <- if(cases) matrix(0, n, m) else NULL
    slack <- numeric(length(boundary.components))
    names(slack) <- names(boundary.rho)
    working <- list()
    for(name in unique(vapply(smooths, `[[`, character(1L),
        "parameter")))
      working[[name]] <- .update(par = par, y = y, eta = ee[[name]],
        family = family, which = name)
    for(z in smooths) {
      ew <- working[[z$parameter]]
      h <- as.numeric(ew$weights)
      if(length(h) != n || any(!is.finite(h)) || any(h <= 0))
        stop("invalid final RS working weights for sandwich covariance.",
          call. = FALSE)
      w <- case.weights * h
      response <- as.numeric(ew$eta) - ee[[z$parameter]] +
        drop(z$X %*% theta[z$index])
      A <- crossprod(z$X, z$X * w)
      lambda <- z$lambda
      lambda[vapply(z$ids, function(i) components[[i]]$penalty,
        integer(1L))] <- exp(r[z$ids])
      Sl <- matrix(0, ncol(z$X), ncol(z$X))
      for(k in seq_along(z$S)) Sl <- Sl + lambda[k] * z$S[[k]]
      L <- tryCatch(chol(A + Sl), error = function(e) NULL)
      if(is.null(L))
        stop("local RS working information is singular.", call. = FALSE)
      Q <- chol2inv(L)
      b <- drop(Q %*% crossprod(z$X, w * response))
      edf <- sum(A * Q)
      residual <- response - drop(z$X %*% b)
      rss <- sum(w * residual^2)
      df <- sum(w > 0) - edf
      if(!is.finite(rss) || rss <= 0 || !is.finite(df) || df <= 0)
        stop("local RS residual scale is not regular.", call. = FALSE)
      sig2 <- rss / df
      Ps <- NULL
      if(length(z$S) > 1L) {
        St <- matrix(0, z$basis$rank, z$basis$rank)
        for(k in seq_along(z$S))
          St <- St + lambda[k] * z$basis$penalties[[k]]
        Ls <- tryCatch(chol(St), error = function(e) NULL)
        if(is.null(Ls))
          stop("local RS penalty range is singular.", call. = FALSE)
        Ps <- chol2inv(Ls)
      }
      if(cases) {
        QX <- Q %*% t(z$X)
        dEDF <- w * colSums(QX * (Sl %*% QX))
        dRSS <- w * residual^2 - 2 * w * residual *
          drop(crossprod(Sl %*% b, QX))
        dscale <- dRSS / rss + dEDF / df
      }
      for(j in seq_along(z$boundary)) {
        k <- z$boundary[j]
        quad <- drop(crossprod(b, z$S[[k]] %*% b))
        if(length(z$S) == 1L) {
          tau2 <- quad / (edf - z$null)
          if(!is.finite(tau2))
            stop("RS local REML boundary update is not finite.",
              call. = FALSE)
          raw <- sig2 / max(tau2, 1e-7)
        } else {
          num <- sum(Ps * z$basis$penalties[[k]]) -
            sum(Q * z$S[[k]])
          if(quad <= 0 && num <= 0) raw <- lambda[k]
          else if(quad <= 0) raw <- Inf
          else if(num <= 0) raw <- 0
          else raw <- lambda[k] * sig2 * num / quad
        }
        slack[z$boundary.ids[j]] <- if(lambda[k] == 1e7)
          log(raw / lambda[k]) else log(lambda[k] / raw)
      }
      for(i in z$ids) {
        k <- components[[i]]$penalty
        num <- if(length(z$S) == 1L) {
          (edf - z$null) / lambda[k]
        } else {
          sum(Ps * z$basis$penalties[[k]]) - sum(Q * z$S[[k]])
        }
        quad <- drop(crossprod(b, z$S[[k]] %*% b))
        if(length(z$S) == 1L) {
          tau2 <- quad / (edf - z$null)
          if(!is.finite(tau2) || !is.finite(num) || num <= 0)
            stop("local RS smoothing equation is nonregular.",
              call. = FALSE)
          if(abs(tau2 / 1e-7 - 1) <= 1e-4)
            stop("local RS smoothing equation is at the variance-floor transition.",
              call. = FALSE)
          if(tau2 < 1e-7) {
            out[i] <- log(sig2) - log(1e-7) - r[i]
            if(cases) C[, i] <- dscale
            next
          }
        }
        if(!is.finite(num) || num <= 0 || !is.finite(quad) ||
            quad <= 0)
          stop("local RS smoothing equation is nonregular.", call. = FALSE)
        out[i] <- log(sig2) + log(num) - log(quad)
        if(cases) {
          dnum <- w * colSums(QX * (z$S[[k]] %*% QX))
          dquad <- 2 * w * residual *
            drop(crossprod(z$S[[k]] %*% b, QX))
          C[, i] <- dscale + dnum / num - dquad / quad
        }
      }
    }
    list(equations = out, cases = C, boundary.slack = slack)
  }

  base <- equations(beta, rho, cases = TRUE)
  if(length(boundary.rho) &&
      (any(is.na(base$boundary.slack)) ||
        any(abs(base$boundary.slack) <= 1e-4)))
    stop("RS smoothing parameter has an indeterminate local REML limit.",
      call. = FALSE)
  if(length(boundary.rho) && any(base$boundary.slack < -1e-4))
    stop("RS fitted smoothing parameter is inconsistent with the final local REML update.",
      call. = FALSE)
  if(m && max(abs(base$equations)) > 0.01)
    stop("RS smoothness equations are not sufficiently close to their fixed point.",
      call. = FALSE)
  p <- length(active)
  q <- p + m
  A <- matrix(0, q, q)
  if(p) A[seq_len(p), seq_len(p)] <- crossprod(R)
  if(m) {
    A[seq_len(p), p + seq_len(m)] <- state$u
    for(j in active) {
      h <- 1e-5 * max(1, abs(beta[j]))
      theta <- beta
      theta[j] <- beta[j] + h
      plus <- equations(theta, rho)$equations
      theta[j] <- beta[j] - h
      minus <- equations(theta, rho)$equations
      A[p + seq_len(m), lookup[j]] <- -(plus - minus) / (2 * h)
    }
    for(j in seq_len(m)) {
      h <- 1e-4
      r <- rho
      r[j] <- rho[j] + h
      plus <- equations(beta, r)$equations
      r[j] <- rho[j] - h
      minus <- equations(beta, r)$equations
      A[p + seq_len(m), p + j] <- -(plus - minus) / (2 * h)
    }
  }
  dimnames(A) <- list(c(names(beta)[active], cn),
    c(names(beta)[active], cn))
  scale <- sqrt(abs(diag(A)))
  condition <- if(q && all(is.finite(scale)) && all(scale > 0))
    rcond(A / outer(scale, scale)) else if(q) 0 else 1
  if(!is.finite(condition) || condition < 1e-10)
    stop("RS sandwich estimating equations are not identifiable.",
      call. = FALSE)

  B <- matrix(0, q, q, dimnames = dimnames(A))
  par <- family$map2par(eta)
  score <- lapply(names(blocks), function(name) {
    s <- as.numeric(family$score[[name]](par = par, y = y))
    if(length(s) != n || any(!is.finite(s)))
      stop("invalid linked score for RS sandwich covariance.",
        call. = FALSE)
    prior * s
  })
  coefficient.residual <- NULL
  if(details) {
    coefficient.residual <- numeric(length(beta))
    for(a in seq_along(blocks))
      coefficient.residual[index[[a]]] <-
        drop(crossprod(X[[a]], score[[a]]))
    coefficient.residual <- coefficient.residual - drop(P %*% beta)
    names(coefficient.residual) <- names(beta)
  }
  for(a in seq_along(blocks)) {
    ia <- index[[a]]
    ka <- which(ia %in% active)
    if(!length(ka)) next
    ia <- lookup[ia[ka]]
    Xa <- X[[a]][, ka, drop = FALSE]
    for(b in a:length(blocks)) {
      ib <- index[[b]]
      kb <- which(ib %in% active)
      if(!length(kb)) next
      ib <- lookup[ib[kb]]
      block <- crossprod(Xa, X[[b]][, kb, drop = FALSE] *
        (score[[a]] * score[[b]]))
      B[ia, ib] <- block
      if(a != b) B[ib, ia] <- t(block)
    }
    if(m) {
      block <- crossprod(Xa, base$cases * score[[a]])
      B[ia, p + seq_len(m)] <- block
      B[p + seq_len(m), ia] <- t(block)
    }
  }
  if(m) B[p + seq_len(m), p + seq_len(m)] <-
    crossprod(base$cases)
  inverse <- if(q) solve(A / outer(scale, scale)) /
    outer(scale, scale) else A
  W <- inverse %*% B %*% t(inverse)
  W <- 0.5 * (W + t(W))
  dimnames(W) <- dimnames(A)
  answer <- V * 0
  if(p) answer[active, active] <- W[seq_len(p), seq_len(p), drop = FALSE]
  if(!details) return(answer)
  Vr <- W[p + seq_len(m), p + seq_len(m), drop = FALSE]
  cross <- matrix(0, length(beta), m,
    dimnames = list(names(beta), cn))
  if(m && p) cross[active, ] <- W[seq_len(p), p + seq_len(m),
    drop = FALSE]
  joint <- rbind(cbind(answer, cross), cbind(t(cross), Vr))
  dimnames(joint) <- list(c(names(beta), cn), c(names(beta), cn))
  list(covariance = answer, smoothing = Vr, cross = cross,
    joint = joint, bread = A, meat = B, rho = rho,
    boundary.rho = boundary.rho,
    boundary.slack = base$boundary.slack,
    smoothing.names = cn, rank = q, condition = condition,
    residual = base$equations,
    coefficient.residual = coefficient.residual[active],
    case.derivative = base$cases,
    equations = equations)
}
