## Joint REML estimation using RS or CG for fixed smoothing parameters.
.JR <- function(x, y, specials, family, offsets, weights, start, xterms, sterms,
  control)
{
  ## Only ordinary quadratic smooth penalties have a joint REML criterion.
  if(isTRUE(control$ridge) || isTRUE(control$termselect))
    stop("JR does not support selected ridge or adaptive term penalties.")
  inner <- control$jr.inner
  if(is.null(inner)) inner <- RS
  if(!identical(inner, RS) && !identical(inner, CG))
    stop("jr.inner must be RS or CG.")
  if(!is.null(control$jr.gradient) &&
      (!is.character(control$jr.gradient) ||
        length(control$jr.gradient) != 1L ||
        is.na(control$jr.gradient) ||
        !control$jr.gradient %in% c("hybrid", "implicit", "finite")))
    stop("jr.gradient must be 'hybrid', 'implicit', or 'finite'.")

  for(j in names(sterms)) for(k in sterms[[j]]) {
    sk <- specials[[k]]
    if(!inherits(sk, "mgcv.smooth") || is.null(sk$X) ||
        is.null(sk$S))
      stop(".JR requires coefficient-linear mgcv smooths with fixed quadratic penalties.")
    if(isTRUE(sk$fixed) || !is.null(sk$sp)) next
    criterion <- sk$control$criterion
    if(is.null(criterion)) criterion <- sk$control$method
    if(is.null(criterion)) criterion <- control$criterion
    if(!is.null(criterion) && tolower(criterion) != "reml")
      stop(".JR requires criterion = 'reml' for every estimated smooth.")
  }

  inner.control <- control
  inner.control$trace <- FALSE
  inner.control$flush <- FALSE
  inner.control$light <- FALSE
  inner.control$.jr.sp <- NULL
  inner.control$.jr.initial.fit <- NULL
  if(is.null(control$rs.cache)) inner.control$rs.cache <- FALSE
  inner.control$eps <- if(is.null(control$jr.eps)) 1e-09 else control$jr.eps
  inner.control$maxit <- if(is.null(control$jr.maxit)) 100L else control$jr.maxit

  ## coef() abbreviates names for one-parameter models by default.
  if(inherits(start, "coef.gamlss2") && length(family$names) == 1L) {
    j <- family$names[1L]
    nm <- names(start)
    if(length(nm) && !any(startsWith(nm, paste0(j, ".")))) {
      linear <- nm %in% xterms[[j]]
      nm[linear] <- paste0(j, ".p.", nm[linear])
      for(k in sterms[[j]]) {
        smooth <- startsWith(nm, paste0(k, "."))
        nm[smooth] <- paste0(j, ".s.", nm[smooth])
      }
      names(start) <- nm
    }
  }

  ## Reuse a complete fit when supplied. Otherwise use local REML to
  ## initialize the outer search.
  fit0 <- control$.jr.initial.fit
  initial <- if(is.null(fit0)) "local" else "fitted"
  if(is.null(fit0)) {
    ## A coefficient start also needs lambdas to avoid the local REML fit.
    if(inherits(start, "coef.gamlss2")) {
      supplied <- unlist(start)
      sp <- setNames(vector("list", length(family$names)), family$names)
      nsp <- nfound <- 0L
      for(j in names(sterms)) for(k in sterms[[j]]) {
        sk <- specials[[k]]
        if(isTRUE(sk$fixed) || !is.null(sk$sp)) next
        tags <- paste0(j, ".s.", k, ".lambda", seq_along(sk$S))
        nsp <- nsp + length(tags)
        if(all(tags %in% names(supplied)) &&
            all(is.finite(supplied[tags])) &&
            all(supplied[tags] > 0)) {
          sp[[j]][[k]] <- as.numeric(supplied[tags])
          nfound <- nfound + length(tags)
        }
      }
      if(nsp > 0L && nfound == nsp) {
        inner.control$.jr.sp <- sp
        initial <- "fixed"
      }
    }
    fit0 <- inner(x, y, specials, family, offsets, weights, start, xterms,
      sterms, inner.control)
  }
  fit0$family <- family
  if(!isTRUE(fit0$converged))
    stop("the initial RS/CG fit did not converge for .JR.")
  unavailable <- function(reason) {
    if(!is.null(fit0[["jr", exact = TRUE]]) &&
        isTRUE(fit0$jr$converged)) {
      fit0$reml.failure <- list(reason = reason, optimizer = ".JR")
      if(isTRUE(control$trace))
        cat("GAMLSS-.JR: joint REML unavailable; previous JR fit retained (",
          reason, ").\n", sep = "")
      return(fit0)
    }
    base <- fit0
    if(identical(initial, "local")) {
      base.control <- control
      base.control$optimizer <- inner
      base.control$trace <- FALSE
      base.control$flush <- FALSE
      base.control$.jr.initial.fit <- NULL
      base.control$.jr.sp <- NULL
      base <- inner(x, y, specials, family, offsets, weights, start,
        xterms, sterms, base.control)
      base$family <- family
    }
    base$jr <- NULL
    base$reml.failure <- list(reason = reason, optimizer = ".JR")
    base$control <- control
    base$control$optimizer <- inner
    base$control$.jr.initial.fit <- NULL
    if(isTRUE(control$trace))
      cat("GAMLSS-.JR: joint REML unavailable; RS/CG fit retained (",
        reason, ").\n", sep = "")
    base
  }

  blocks <- list()
  rho0 <- numeric(0L)
  for(j in names(sterms)) for(k in sterms[[j]]) {
    sk <- specials[[k]]
    sf <- fit0$fitted.specials[[j]][[k]]
    if(!inherits(sk, "mgcv.smooth") || is.null(sk$X) ||
        is.null(sk$S) || is.null(sf$coefficients) ||
        !is.null(sf$penalty))
      stop(".JR requires coefficient-linear mgcv smooths with fixed quadratic penalties.")
    S <- lapply(sk$S, function(z) {
      z <- as.matrix(z)
      0.5 * (z + t(z))
    })
    rank <- if(length(sk$null.space.dim) == 1L)
      ncol(sk$X) - sk$null.space.dim else NULL
    basis <- if(length(S)) smooth.construct_reml(S, rank) else NULL
    base.logdet <- if(!is.null(basis) && length(S) == 1L)
      2 * sum(log(diag(chol(basis$penalties[[1L]])))) else NULL
    lambda <- as.numeric(sf$lambdas)
    if(length(lambda) != length(S) || any(!is.finite(lambda)) ||
        any(lambda <= 0))
      stop(".JR requires positive fitted smoothing parameters.")
    fixed <- isTRUE(sk$fixed) || !is.null(sk$sp)
    ids <- if(fixed) integer(0L) else
      seq.int(length(rho0) + 1L, length(rho0) + length(S))
    if(length(ids)) {
      rho0 <- c(rho0, log(lambda))
      names(rho0)[ids] <- paste0(j, ".s.", k, ".lambda", seq_along(S))
    }
    blocks[[length(blocks) + 1L]] <- list(parameter = j, term = k,
      S = S, basis = basis, base.logdet = base.logdet,
      fixed = fixed, ids = ids, lambda = lambda)
  }

  if(!length(rho0)) {
    fit0$jr <- list(converged = TRUE, initial = initial, rho = rho0,
      covariance = matrix(0, 0L, 0L))
    fit0$control <- control
    fit0$control$.jr.initial.fit <- NULL
    return(fit0)
  }

  ## Reuse the model setup and warm starts. Only the fixed smoothing
  ## parameters and coefficients change between outer evaluations.
  state <- new.env(parent = emptyenv())
  state$best <- fit0
  state$value <- Inf
  state$last.rho <- NULL
  state$last.fit <- NULL
  state$last.V <- NULL
  state$last.error <- NULL
  state$active <- NULL
  state$last.value <- Inf
  state$evaluations <- 0L
  state$gradient.design <- NULL
  state$gradient.fallbacks <- 0L

  attach_setup <- function(fit) {
    fit$family <- family
    fit$x <- x
    fit$y <- y
    fit$specials <- specials
    fit$xterms <- xterms
    fit$sterms <- sterms
    fit$weights <- weights
    fit$offsets <- offsets
    fit
  }

  criterion <- function(rho) {
    if(!is.null(state$last.rho) && identical(rho, state$last.rho))
      return(state$last.value)
    state$evaluations <- state$evaluations + 1L
    state$last.error <- NULL
    sp <- setNames(vector("list", length(family$names)), family$names)
    initial <- coef(state$best, full = TRUE, lambdas = TRUE,
      dropall = FALSE)
    for(z in blocks) if(length(z$ids)) {
      lambda <- exp(rho[z$ids])
      sp[[z$parameter]][[z$term]] <- lambda
      tags <- paste0(z$parameter, ".s.", z$term, ".lambda",
        seq_along(lambda))
      initial[tags] <- lambda
    }
    trial.control <- inner.control
    trial.control$.jr.sp <- sp
    out <- tryCatch({
      fit <- inner(x, y, specials, family, offsets, weights, initial,
        xterms, sterms, trial.control)
      if(!isTRUE(fit$converged))
        stop("inner fit did not converge")
      fit <- attach_setup(fit)
      V <- vcov.gamlss2(fit, full = TRUE, method = "joint")
      if(anyNA(diag(V)))
        stop("joint information is singular")
      active <- which(diag(V) > 0)
      if(!length(active) || any(!is.finite(V[active, active])))
        stop("joint information is singular")
      active.names <- rownames(V)[active]
      if(!is.null(state$active) && !identical(active.names, state$active))
        stop("the estimable coefficient set changed")
      V <- V[active, active, drop = FALSE]
      R <- chol(V)
      logdetK <- -2 * sum(log(diag(R)))
      penalty <- logdetP <- 0
      for(z in blocks) {
        b <- as.numeric(fit$fitted.specials[[z$parameter]][[z$term]]$coefficients)
        lambda <- as.numeric(fit$fitted.specials[[z$parameter]][[z$term]]$lambdas)
        P <- matrix(0, length(b), length(b))
        for(i in seq_along(z$S)) P <- P + lambda[i] * z$S[[i]]
        penalty <- penalty + drop(crossprod(b, P %*% b))
        if(!is.null(z$basis)) {
          if(length(z$S) == 1L) {
            logdetP <- logdetP + z$base.logdet +
              z$basis$rank * log(lambda[1L])
          } else {
            Pr <- matrix(0, z$basis$rank, z$basis$rank)
            for(i in seq_along(z$S))
              Pr <- Pr + lambda[i] * z$basis$penalties[[i]]
            logdetP <- logdetP + 2 * sum(log(diag(chol(Pr))))
          }
        }
      }
      value <- -fit$logLik + 0.5 * (penalty + logdetK - logdetP)
      if(!is.finite(value)) stop("joint REML is not finite")
      state$active <- active.names
      list(value = value, fit = fit, V = V)
    }, error = function(e) {
      state$last.error <- conditionMessage(e)
      NULL
    })
    value <- if(is.null(out)) Inf else out$value
    if(is.finite(value) && value < state$value) {
      state$value <- value
      state$best <- out$fit
    }
    state$last.rho <- rho
    state$last.fit <- if(is.null(out)) NULL else out$fit
    state$last.V <- if(is.null(out)) NULL else out$V
    state$last.value <- value
    value
  }

  lower <- rep(log(1e-10), length(rho0))
  upper <- rep(log(1e+10), length(rho0))
  gradient.mode <- control$jr.gradient
  if(is.null(gradient.mode))
    gradient.mode <- if(length(rho0) < 3L) "finite" else "hybrid"
  if(!is.finite(criterion(rho0)))
    return(unavailable(paste0(
      "initial joint REML evaluation failed: ", state$last.error)))

  ## The coefficient derivative follows from the penalized score equation.
  ## Only the likelihood-information determinant needs a directional
  ## difference; this avoids refitting coefficients for each gradient.
  gradient_design <- function(fit) {
    p <- length(state$active)
    n <- if(is.null(dim(y))) length(y) else nrow(y)
    D <- setNames(lapply(family$names, function(j)
      matrix(0, n, p)), family$names)
    indices <- vector("list", length(blocks))
    co <- coef(fit, full = TRUE, dropall = FALSE)
    for(j in family$names) {
      b <- fit$coefficients[[j]]
      if(length(b)) {
        ii <- match(paste0(j, ".p.", names(b)), state$active)
        keep <- which(!is.na(ii))
        if(length(keep))
          D[[j]][, ii[keep]] <- x[, names(b)[keep], drop = FALSE]
      }
    }
    for(a in seq_along(blocks)) {
      z <- blocks[[a]]
      sk <- specials[[z$term]]
      X <- sk$X
      if(nrow(X) != n && !is.null(sk$binning))
        X <- X[sk$binning$match.index, , drop = FALSE]
      prefix <- paste0(z$parameter, ".s.", z$term, ".")
      tags <- names(co)[startsWith(names(co), prefix)]
      if(nrow(X) != n || length(tags) != ncol(X))
        stop("cannot construct the joint REML gradient design")
      ii <- match(tags, state$active)
      indices[[a]] <- ii
      keep <- which(!is.na(ii))
      if(length(keep))
        D[[z$parameter]][, ii[keep]] <- X[, keep, drop = FALSE]
    }
    list(D = D, indices = indices)
  }

  implicit_gradient <- function(rho) {
    criterion(rho)
    fit <- state$last.fit
    V <- state$last.V
    if(is.null(fit) || is.null(V))
      stop("joint REML information is unavailable")
    if(is.null(state$gradient.design))
      state$gradient.design <- gradient_design(fit)
    D <- state$gradient.design$D
    eta <- as.data.frame(fit$fitted.values)
    beta <- as.numeric(coef(fit, full = TRUE, dropall = FALSE)[state$active])
    g <- numeric(length(rho))
    for(a in seq_along(blocks)) {
      z <- blocks[[a]]
      if(!length(z$ids)) next
      ii <- state$gradient.design$indices[[a]]
      if(anyNA(ii))
        stop("joint REML smooth coefficients are not estimable")
      lambda <- exp(rho[z$ids])
      Pr <- NULL
      if(!is.null(z$basis) && length(z$S) > 1L) {
        Pr <- matrix(0, z$basis$rank, z$basis$rank)
        for(k in seq_along(z$S))
          Pr <- Pr + lambda[k] * z$basis$penalties[[k]]
        Pr <- chol2inv(chol(Pr))
      }
      for(k in seq_along(z$ids)) {
        Pi <- lambda[k] * z$S[[k]]
        u <- numeric(length(beta))
        u[ii] <- drop(Pi %*% beta[ii])
        db <- -drop(V %*% u)
        deta <- lapply(D, function(X) drop(X %*% db))
        scale <- max(abs(unlist(deta, use.names = FALSE)))
        h <- if(scale > 0) min(0.05, 0.01 / scale) else 0.05
        logdet <- function(sign) {
          trial <- fit
          trial$fitted.values <- eta
          for(j in names(D))
            trial$fitted.values[[j]] <- eta[[j]] + sign * h * deta[[j]]
          W <- vcov.gamlss2(trial, full = TRUE, method = "joint")
          W <- W[state$active, state$active, drop = FALSE]
          -2 * sum(log(diag(chol(W))))
        }
        penalty.trace <- if(is.null(z$basis)) 0 else
          if(length(z$S) == 1L) z$basis$rank else
            sum(Pr * (lambda[k] * z$basis$penalties[[k]]))
        g[z$ids[k]] <- 0.5 * (sum(beta[ii] * u[ii]) +
          sum(V[ii, ii, drop = FALSE] * Pi) - penalty.trace +
          (logdet(1) - logdet(-1)) / (2 * h))
      }
    }
    if(any(!is.finite(g))) stop("joint REML gradient is not finite")
    g
  }

  finite_gradient <- function(rho) {
    h <- 0.01
    g <- numeric(length(rho))
    for(i in seq_along(rho)) {
      r <- rho
      r[i] <- min(upper[i], rho[i] + h)
      fp <- criterion(r)
      r[i] <- max(lower[i], rho[i] - h)
      fm <- criterion(r)
      g[i] <- (fp - fm) / (min(upper[i], rho[i] + h) -
        max(lower[i], rho[i] - h))
    }
    g
  }
  outer_gradient <- if(identical(gradient.mode, "finite")) {
    finite_gradient
  } else {
    function(rho) {
      g <- tryCatch(implicit_gradient(rho), error = function(e) NULL)
      if(!is.null(g)) return(g)
      state$gradient.fallbacks <- state$gradient.fallbacks + 1L
      finite_gradient(rho)
    }
  }
  opt <- tryCatch(nlminb(rho0, criterion, gradient = outer_gradient,
    lower = lower, upper = upper,
    control = list(iter.max = if(is.null(control$jr.outer.maxit)) 30L else
      control$jr.outer.maxit, eval.max = if(is.null(control$jr.outer.maxeval))
      150L else control$jr.outer.maxeval, rel.tol = 1e-08)),
    error = function(e) e)
  if(inherits(opt, "error"))
    return(unavailable(conditionMessage(opt)))
  if(identical(gradient.mode, "hybrid") && opt$convergence == 0L) {
    first <- opt
    polished <- tryCatch(nlminb(opt$par, criterion,
      gradient = finite_gradient,
      lower = lower, upper = upper,
      control = list(iter.max = 5L, eval.max = 40L,
        rel.tol = 1e-08)), error = function(e) NULL)
    if(!is.null(polished) && is.finite(polished$objective) &&
        polished$objective <= first$objective + 1e-10) {
      opt <- polished
      if(opt$convergence != 0L)
        opt$convergence <- first$convergence
      opt$iterations <- first$iterations + polished$iterations
      opt$evaluations <- first$evaluations + polished$evaluations
      opt$initial <- first
    }
  }
  rho <- opt$par
  value <- criterion(rho)
  fit <- state$last.fit
  if(!is.finite(value) || is.null(fit))
    return(unavailable("optimized joint REML evaluation failed"))

  ## A small deviance change alone does not establish a coefficient fixed
  ## point. Check one more sweep at the final smoothing parameters.
  check.sp <- setNames(vector("list", length(family$names)), family$names)
  for(z in blocks) if(length(z$ids))
    check.sp[[z$parameter]][[z$term]] <- exp(rho[z$ids])
  check.control <- inner.control
  check.control$.jr.sp <- check.sp
  check.control$maxit <- 1L
  check <- inner(x, y, specials, family, offsets, weights,
    coef(fit, full = TRUE, lambdas = TRUE, dropall = FALSE),
    xterms, sterms, check.control)
  check$family <- family
  b0 <- as.numeric(coef(fit, full = TRUE, dropall = FALSE))
  b1 <- as.numeric(coef(check, full = TRUE, dropall = FALSE))
  residual <- max(abs(b1 - b0)) / max(1, abs(b0))

  ## Differentiate the profiled criterion, including coefficient refits.
  ## Smoothing parameters with vanishing curvature are treated as being
  ## at their effective boundary.
  m <- length(rho)
  h <- rep(0.01, m)
  gradient <- rep(NA_real_, m)
  H <- matrix(NA_real_, m, m)
  dC <- vector("list", m)
  covariance_factor <- function(V) {
    R <- chol(chol2inv(chol(V)))
    backsolve(R, diag(nrow(V)))
  }
  interior <- all(rho - h > lower & rho + h < upper)
  if(interior) {
    plus <- minus <- numeric(m)
    for(i in seq_len(m)) {
      r <- rho
      r[i] <- r[i] + h[i]
      plus[i] <- criterion(r)
      Vp <- state$last.V
      r[i] <- rho[i] - h[i]
      minus[i] <- criterion(r)
      Vm <- state$last.V
      gradient[i] <- (plus[i] - minus[i]) / (2 * h[i])
      H[i, i] <- (plus[i] - 2 * value + minus[i]) / h[i]^2
      if(!is.null(Vp) && !is.null(Vm)) {
        Cp <- tryCatch(covariance_factor(Vp), error = function(e) NULL)
        Cm <- tryCatch(covariance_factor(Vm), error = function(e) NULL)
        if(!is.null(Cp) && !is.null(Cm))
          dC[[i]] <- (Cp - Cm) / (2 * h[i])
      }
    }
    if(m > 1L) for(i in seq_len(m - 1L)) for(j in (i + 1L):m) {
      v <- numeric(4L)
      signs <- rbind(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
      for(k in seq_len(4L)) {
        r <- rho
        r[i] <- r[i] + signs[k, 1L] * h[i]
        r[j] <- r[j] + signs[k, 2L] * h[j]
        v[k] <- criterion(r)
      }
      H[i, j] <- H[j, i] <- (v[1L] - v[2L] - v[3L] + v[4L]) /
        (4 * h[i] * h[j])
    }
  }
  covariance <- NULL
  factor.root <- NULL
  curvature.rank <- 0L
  curvature.cutoff <- NA_real_
  curvature.values <- rep(NA_real_, m)
  stationary <- FALSE
  if(all(is.finite(H)) && all(is.finite(gradient))) {
    ev <- eigen(0.5 * (H + t(H)), symmetric = TRUE)
    curvature.values <- ev$values
    curvature.cutoff <- max(1e-05 * max(abs(ev$values)),
      100 * .Machine$double.eps * max(1, abs(value)) / min(h)^2)
    positive <- ev$values > curvature.cutoff
    curvature.rank <- sum(positive)
    if(any(positive) && !any(ev$values < -curvature.cutoff)) {
      Q <- ev$vectors[, positive, drop = FALSE]
      Q <- sweep(Q, 2L, sqrt(ev$values[positive]), "/")
      covariance <- tcrossprod(Q)
      step <- drop(covariance %*% gradient)
      flat.gradient <- gradient -
        drop(ev$vectors[, positive, drop = FALSE] %*%
          crossprod(ev$vectors[, positive, drop = FALSE], gradient))
      stationary <- max(abs(step)) < 0.05 &&
        sum(gradient * step) < 0.04 &&
        max(abs(flat.gradient)) < 1e-05 &&
        all(vapply(dC, function(x) !is.null(x) && all(is.finite(x)),
          logical(1L)))
      ## As in mgcv, bound covariance-factor variation in flat directions.
      ## The mean correction uses only the identified curvature subspace.
      factor.variance <- rep(if(length(family$names) > 1L) 50 else 10, m)
      factor.variance[positive] <- 1 / ev$values[positive]
      factor.root <- t(sweep(ev$vectors, 2L,
        sqrt(factor.variance), "*"))
    }
  }
  if(!stationary) {
    covariance <- NULL
    factor.root <- NULL
    dC <- NULL
  }
  names(rho) <- names(rho0)
  names(gradient) <- names(rho0)
  if(!is.null(dC)) names(dC) <- names(rho0)
  dimnames(H) <- list(names(rho0), names(rho0))
  if(!is.null(covariance))
    dimnames(covariance) <- dimnames(H)
  if(!is.null(factor.root))
    colnames(factor.root) <- names(rho0)
  fit$jr <- list(initial = initial, rho = rho, criterion = value,
    gradient = gradient,
    hessian = H, covariance = covariance, factor.root = factor.root,
    curvature.rank = curvature.rank, curvature.cutoff = curvature.cutoff,
    curvature.values = curvature.values, factor.derivative = dC,
    fixed.point.residual = residual,
    converged = opt$convergence == 0L && stationary &&
      is.finite(residual) && residual < 1e-05,
    evaluations = state$evaluations, gradient.mode = gradient.mode,
    gradient.fallbacks = state$gradient.fallbacks, outer = opt)
  fit$control <- control
  fit$control$.jr.initial.fit <- NULL
  fit$converged <- isTRUE(fit$converged) && fit$jr$converged
  if(!fit$jr$converged)
    return(unavailable(
      "joint REML curvature or coefficient fixed point could not be verified"))
  if(isTRUE(control$trace))
    cat("GAMLSS-.JR: joint REML = ", signif(value, 8L),
      " evaluations = ", state$evaluations,
      " converged = ", fit$jr$converged, "\n", sep = "")
  fit
}


## Continue joint REML estimation from a complete fitted model.
.jr <- function(object, trace = TRUE, ...)
{
  if(!inherits(object, "gamlss2") || inherits(object, "bamlss2"))
    stop(".jr requires a fitted gamlss2 model.")
  if(is.null(object$x) || is.null(object$y) ||
      is.null(object$specials) || is.null(object$xterms) ||
      is.null(object$sterms))
    stop(".jr requires a fit retaining x, y, and the smooth setup.")
  if(any(vapply(object$specials, function(z)
      inherits(z, "mgcv.smooth") && is.null(z$X), logical(1L))))
    stop(".jr requires the stored smooth design matrices.")

  control <- object$control
  control$optimizer <- .JR
  control$trace <- trace
  extra <- list(...)
  if(length(extra)) {
    if(is.null(names(extra)) || any(!nzchar(names(extra))))
      stop(".jr control arguments must be named.")
    for(k in names(extra)) control[[k]] <- extra[[k]]
  }
  control$.jr.initial.fit <- object

  tstart <- proc.time()
  fit <- .JR(x = object$x, y = object$y, specials = object$specials,
    family = object$family, offsets = object$offsets,
    weights = object$weights,
    start = coef(object, full = TRUE, lambdas = TRUE, dropall = FALSE),
    xterms = object$xterms, sterms = object$sterms, control = control)
  elapsed <- as.numeric((proc.time() - tstart)["elapsed"])
  if(is.null(fit[["jr", exact = TRUE]])) {
    object$reml.failure <- fit$reml.failure
    object$reml.failure$elapsed <- elapsed
    return(object)
  }

  ## Keep the original model call and setup for prediction and plotting.
  object$reml.failure <- NULL
  for(k in names(fit)) object[[k]] <- fit[[k]]
  object$df <- get_df(object)
  object$elapsed <- object$elapsed + elapsed
  object$jr$elapsed <- elapsed
  object$jr.call <- match.call()
  if(!is.null(object$call))
    object$call$optimizer <- quote(.JR)
  if(!is.null(object$model))
    object$results <- results(object, data = object$model, interval = "local")
  object
}


## Joint REML with a cached design and a direct coefficient solver.
JR <- function(x, y, specials, family, offsets, weights, start, xterms,
  sterms, control)
{
  kind <- family$family
  links <- unname(unlist(family$links))
  analytic <- (identical(kind, "NO") && identical(links,
    c("identity", "log"))) ||
    (identical(kind, "GA") && identical(links, c("log", "log"))) ||
    (identical(kind, "PO") && identical(links, "log"))
  if(isTRUE(control$ridge) || isTRUE(control$termselect) ||
      any(vapply(control$fixed, isTRUE, logical(1L))))
    stop("JR requires ordinary estimable coefficients and quadratic smooth penalties.")
  if(!is.function(family$map2par) || !is.function(family$pdf))
    stop("JR requires a family density and predictor-to-parameter map.")

  inner <- control$jr.inner
  if(is.null(inner)) inner <- RS
  if(!identical(inner, RS) && !identical(inner, CG))
    stop("jr.inner must be RS or CG.")
  for(j in names(sterms)) for(k in sterms[[j]]) {
    sk <- specials[[k]]
    if(!inherits(sk, "mgcv.smooth") || is.null(sk$X) ||
        is.null(sk$S))
      stop("JR requires coefficient-linear mgcv smooths with quadratic penalties.")
    if(isTRUE(sk$fixed) || !is.null(sk$sp)) next
    criterion <- sk$control$criterion
    if(is.null(criterion)) criterion <- sk$control$method
    if(is.null(criterion)) criterion <- control$criterion
    if(!is.null(criterion) && tolower(criterion) != "reml")
      stop("JR requires criterion = 'reml' for every estimated smooth.")
  }

  inner.control <- control
  inner.control$trace <- FALSE
  inner.control$flush <- FALSE
  inner.control$light <- FALSE
  inner.control$.jr.sp <- NULL
  inner.control$.jr.initial.fit <- NULL
  if(is.null(control$rs.cache)) inner.control$rs.cache <- FALSE
  inner.control$eps <- if(is.null(control$jr.eps)) 1e-09 else control$jr.eps
  inner.control$maxit <- if(is.null(control$jr.maxit)) 100L else
    control$jr.maxit
  fit0 <- control$.jr.initial.fit
  initial <- if(is.null(fit0)) "local" else "fitted"
  if(is.null(fit0) && inherits(start, "coef.gamlss2")) {
    supplied <- unlist(start)
    np <- family$names
    sp0 <- setNames(vector("list", length(np)), np)
    nsp <- nfound <- 0L
    for(j in names(sterms)) for(k in sterms[[j]]) {
      sk <- specials[[k]]
      if(isTRUE(sk$fixed) || !is.null(sk$sp)) next
      tags <- paste0(j, ".s.", k, ".lambda", seq_along(sk$S))
      nsp <- nsp + length(tags)
      if(all(tags %in% names(supplied)) &&
          all(is.finite(supplied[tags])) &&
          all(supplied[tags] > 0)) {
        sp0[[j]][[k]] <- as.numeric(supplied[tags])
        nfound <- nfound + length(tags)
      }
    }
    if(nsp > 0L && nfound == nsp) {
      inner.control$.jr.sp <- sp0
      initial <- "fixed"
    }
  }
  if(is.null(fit0)) {
    initial.control <- inner.control
    initial.control$eps <- if(is.null(control$eps)) 1e-05 else
      control$eps
    initial.control$maxit <- if(is.null(control$maxit)) 20L else
      control$maxit
    initial.control$rs.cache <- control$rs.cache
    fit0 <- inner(x, y, specials, family, offsets, weights, start, xterms,
      sterms, initial.control)
  }
  fit0$family <- family
  if(!isTRUE(fit0$converged))
    stop("the initial RS/CG fit did not converge for JR.")
  unavailable <- function(reason) {
    if(!is.null(fit0[["jr", exact = TRUE]]) &&
        isTRUE(fit0$jr$converged)) {
      fit0$reml.failure <- list(reason = reason, optimizer = "JR")
      if(isTRUE(control$trace))
        cat("GAMLSS-JR: joint REML unavailable; previous JR fit retained (",
          reason, ").\n", sep = "")
      return(fit0)
    }
    base <- fit0
    if(identical(initial, "local")) {
      base.control <- control
      base.control$optimizer <- inner
      base.control$trace <- FALSE
      base.control$flush <- FALSE
      base.control$.jr.initial.fit <- NULL
      base.control$.jr.sp <- NULL
      base <- inner(x, y, specials, family, offsets, weights, start,
        xterms, sterms, base.control)
      base$family <- family
    }
    base$jr <- NULL
    base$reml.failure <- list(reason = reason, optimizer = "JR")
    base$control <- control
    base$control$optimizer <- inner
    base$control$.jr.initial.fit <- NULL
    if(isTRUE(control$trace))
      cat("GAMLSS-JR: joint REML unavailable; RS/CG fit retained (",
        reason, ").\n", sep = "")
    base
  }

  co <- coef(fit0, full = TRUE, dropall = FALSE)
  beta0 <- as.numeric(co)
  p <- length(beta0)
  n <- if(is.null(dim(y))) length(y) else nrow(y)
  np <- family$names
  if(any(!is.finite(beta0)) ||
      (is.numeric(y) && any(!is.finite(y))))
    stop("JR requires finite initial coefficients and response values.")
  w <- if(is.null(weights)) rep(1, n) else as.numeric(weights)
  if(length(w) != n || any(!is.finite(w)) || any(w < 0))
    stop("invalid observation weights.")
  index <- X <- off <- setNames(vector("list", length(np)), np)
  for(j in np) {
    index[[j]] <- which(startsWith(names(co), paste0(j, ".")))
    X[[j]] <- matrix(0, n, length(index[[j]]),
      dimnames = list(NULL, names(co)[index[[j]]]))
    off[[j]] <- if(is.null(offsets) || !j %in% names(offsets))
      rep(0, n) else as.numeric(offsets[[j]])
    if(length(off[[j]]) != n) stop("invalid model offsets.")
    bj <- fit0$coefficients[[j]]
    if(length(bj)) {
      ii <- match(paste0(j, ".p.", names(bj)), colnames(X[[j]]))
      if(anyNA(ii))
        stop("JR could not match the linear coefficient design.")
      X[[j]][, ii] <- x[, names(bj), drop = FALSE]
    }
  }

  ## Assemble each penalty only once; its scale changes in the outer fit.
  blocks <- list()
  rho0 <- numeric(0L)
  for(j in names(sterms)) for(k in sterms[[j]]) {
    sk <- specials[[k]]
    sf <- fit0$fitted.specials[[j]][[k]]
    Z <- sk$X
    if(nrow(Z) != n && !is.null(sk$binning))
      Z <- Z[sk$binning$match.index, , drop = FALSE]
    tags <- paste0(j, ".s.", k, ".", seq_len(ncol(Z)))
    ii <- match(tags, colnames(X[[j]]))
    bi <- match(tags, names(co))
    if(nrow(Z) != n || anyNA(ii) || anyNA(bi) ||
        !is.null(sf$penalty))
      stop("JR could not match the stored smooth coefficient design.")
    X[[j]][, ii] <- Z
    S <- lapply(sk$S, function(z) {
      z <- as.matrix(z)
      0.5 * (z + t(z))
    })
    lambda <- as.numeric(sf$lambdas)
    if(length(lambda) != length(S) || any(!is.finite(lambda)) ||
        any(lambda <= 0))
      stop("JR requires positive fitted smoothing parameters.")
    rank <- if(length(sk$null.space.dim) == 1L)
      ncol(Z) - sk$null.space.dim else NULL
    basis <- if(length(S)) smooth.construct_reml(S, rank) else NULL
    base.logdet <- if(!is.null(basis) && length(S) == 1L)
      2 * sum(log(diag(chol(basis$penalties[[1L]])))) else NULL
    fixed <- isTRUE(sk$fixed) || !is.null(sk$sp)
    ids <- if(fixed) integer(0L) else
      seq.int(length(rho0) + 1L, length(rho0) + length(S))
    if(length(ids)) {
      rho0 <- c(rho0, log(lambda))
      names(rho0)[ids] <- paste0(j, ".s.", k, ".lambda",
        seq_along(S))
    }
    blocks[[length(blocks) + 1L]] <- list(parameter = j, term = k,
      index = bi, S = S, basis = basis, base.logdet = base.logdet,
      ids = ids, lambda = lambda)
  }
  if(!length(rho0)) {
    fit0$jr <- list(converged = TRUE, initial = initial, rho = rho0,
      covariance = matrix(0, 0L, 0L))
    fit0$control <- control
    fit0$control$.jr.initial.fit <- NULL
    return(fit0)
  }

  numerical.score <- setNames(!vapply(np, function(j)
    is.function(family$score[[j]]), logical(1L)), np)
  linked_score <- function(eta, j, par = NULL) {
    if(!numerical.score[j]) {
      if(is.null(par)) par <- family$map2par(eta)
      return(as.numeric(family$score[[j]](par = par, y = y)))
    }
    step <- 1e-04 * pmax(1, abs(eta[[j]]))
    upper <- lower <- eta
    upper[[j]] <- eta[[j]] + step
    lower[[j]] <- eta[[j]] - step
    fu <- family$pdf(par = family$map2par(upper), y = y, log = TRUE)
    fl <- family$pdf(par = family$map2par(lower), y = y, log = TRUE)
    as.numeric((fu - fl) / (2 * step))
  }

  ## Negative log likelihood, linked-scale score and observed information.
  likelihood <- function(beta, information = TRUE) {
    eta <- lapply(np, function(j)
      drop(X[[j]] %*% beta[index[[j]]]) + off[[j]])
    names(eta) <- np
    if(analytic && kind == "NO") {
      e <- y - eta[[1L]]
      v <- exp(-2 * eta[[2L]])
      q <- e^2 * v
      f <- eta[[2L]] + 0.5 * (q + log(2 * pi))
      g <- list(-e * v, 1 - q)
      if(identical(information, TRUE))
        h <- list(list(v, 2 * e * v), list(2 * e * v, 2 * q))
    } else if(analytic && kind == "GA") {
      a <- exp(-2 * eta[[2L]])
      t <- y * exp(-eta[[1L]])
      b <- digamma(a) - log(a) - 1 + eta[[1L]] + t - log(y)
      f <- lgamma(a) - a * log(a) + a * (eta[[1L]] + t) -
        (a - 1) * log(y)
      g <- list(a * (1 - t), -2 * a * b)
      if(identical(information, TRUE)) {
        cross <- -2 * a * (1 - t)
        h <- list(list(a * t, cross),
          list(cross, 4 * a * b + 4 * a^2 * trigamma(a) - 4 * a))
      }
    } else if(analytic && kind == "PO") {
      mu <- exp(eta[[1L]])
      f <- mu - y * eta[[1L]] + lgamma(y + 1)
      g <- list(mu - y)
      if(identical(information, TRUE)) h <- list(list(mu))
    } else {
      par <- tryCatch(family$map2par(eta), error = function(e) NULL)
      f <- if(is.null(par)) rep(Inf, n) else tryCatch(
        -as.numeric(family$pdf(par = par, y = y, log = TRUE)),
        error = function(e) rep(Inf, n))
      if(length(f) != n) f <- rep(Inf, n)
      if(!identical(information, FALSE) && all(is.finite(f))) {
        g <- lapply(np, function(j) tryCatch(
          -linked_score(eta, j, par = par),
          error = function(e) rep(NA_real_, n)))
        if(any(vapply(g, function(z)
            length(z) != n || any(!is.finite(z)), logical(1L))))
          g <- rep(list(rep(NA_real_, n)), length(np))
        if(identical(information, TRUE)) {
          h <- lapply(np, function(j) vector("list", length(np)))
          for(b in seq_along(np)) {
            step <- 1e-03 * pmax(1, abs(eta[[b]]))
            upper <- lower <- eta
            upper[[b]] <- eta[[b]] + step
            lower[[b]] <- eta[[b]] - step
            pu <- tryCatch(family$map2par(upper),
              error = function(e) NULL)
            pl <- tryCatch(family$map2par(lower),
              error = function(e) NULL)
            for(a in seq_len(b)) {
              su <- if(is.null(pu)) rep(NA_real_, n) else tryCatch(
                linked_score(upper, np[a], par = pu),
                error = function(e) rep(NA_real_, n))
              sl <- if(is.null(pl)) rep(NA_real_, n) else tryCatch(
                linked_score(lower, np[a], par = pl),
                error = function(e) rep(NA_real_, n))
              h[[a]][[b]] <- if(length(su) == n && length(sl) == n)
                -(su - sl) / (2 * step) else rep(NA_real_, n)
            }
          }
        }
      }
    }
    value <- sum(w * f)
    if(identical(information, FALSE)) return(value)
    if(!is.finite(value) || any(!is.finite(unlist(g))) ||
        (identical(information, TRUE) &&
          any(!is.finite(unlist(h)))))
      return(list(value = Inf, score = rep(NA_real_, p),
        K = matrix(NA_real_, p, p)))
    score <- numeric(p)
    for(a in seq_along(np)) {
      ia <- index[[a]]
      score[ia] <- drop(crossprod(X[[a]], w * g[[a]]))
    }
    if(identical(information, "score"))
      return(list(value = value, score = score))
    K <- matrix(0, p, p)
    for(a in seq_along(np)) for(b in a:length(np)) {
      ia <- index[[a]]
      ib <- index[[b]]
      K[ia, ib] <- crossprod(X[[a]], X[[b]] * (w * h[[a]][[b]]))
      if(a != b) K[ib, ia] <- t(K[ia, ib])
    }
    list(value = value, score = score, K = K)
  }

  ## User-modified family callbacks use the general direct solver.
  eta0 <- as.data.frame(lapply(np, function(j)
    drop(X[[j]] %*% beta0[index[[j]]]) + off[[j]]))
  names(eta0) <- np
  check <- tryCatch({
    par <- family$map2par(eta0)
    -sum(family$pdf(par = par, y = y, log = TRUE) * w)
  }, error = function(e) NA_real_)
  if(!is.finite(check))
    stop("JR requires a finite family log likelihood at the initial fit.")
  if(analytic && (!is.finite(likelihood(beta0, FALSE)) ||
      abs(check - likelihood(beta0, FALSE)) >
        1e-07 * max(1, abs(check))))
    analytic <- FALSE
  if(!analytic) {
    eta0 <- as.list(eta0)
    for(j in np) {
      analytic.score <- if(is.function(family$score[[j]])) tryCatch(
        as.numeric(family$score[[j]](
          par = family$map2par(eta0), y = y)),
        error = function(e) rep(NA_real_, n)) else rep(NA_real_, n)
      numerical.score[j] <- TRUE
      density.score <- tryCatch(linked_score(eta0, j),
        error = function(e) rep(NA_real_, n))
      if(length(density.score) != n || any(!is.finite(density.score)))
        stop("JR could not evaluate the family density derivative for parameter '",
          j, "'.")
      if(length(analytic.score) == n && all(is.finite(analytic.score)) &&
          max(abs(analytic.score - density.score) /
            pmax(1, abs(density.score))) < 1e-03 &&
          abs(sum(w * (analytic.score - density.score))) <
            1e-04 * max(1, abs(sum(w * density.score))))
        numerical.score[j] <- FALSE
    }
    if(any(numerical.score)) numerical.score[] <- TRUE
  }

  state <- new.env(parent = emptyenv())
  state$beta <- beta0
  state$last.rho <- NULL
  state$last.value <- Inf
  state$last.beta <- NULL
  state$last.V <- NULL
  state$evaluations <- 0L
  lower <- rep(log(1e-10), length(rho0))
  upper <- rep(log(1e+10), length(rho0))

  profile <- function(rho) {
    if(!is.null(state$last.rho) && identical(rho, state$last.rho))
      return(state$last.value)
    state$evaluations <- state$evaluations + 1L
    P <- matrix(0, p, p)
    logdetP <- 0
    for(z in blocks) {
      lambda <- if(length(z$ids)) exp(rho[z$ids]) else z$lambda
      for(i in seq_along(z$S))
        P[z$index, z$index] <- P[z$index, z$index] +
          lambda[i] * z$S[[i]]
      if(!is.null(z$basis)) {
        if(length(z$S) == 1L) {
          logdetP <- logdetP + z$base.logdet +
            z$basis$rank * log(lambda[1L])
        } else {
          Pr <- matrix(0, z$basis$rank, z$basis$rank)
          for(i in seq_along(z$S))
            Pr <- Pr + lambda[i] * z$basis$penalties[[i]]
          logdetP <- logdetP + 2 * sum(log(diag(chol(Pr))))
        }
      }
    }
    beta <- state$beta
    good <- FALSE
    for(iter in seq_len(30L)) {
      lk <- likelihood(beta)
      if(!is.finite(lk$value) || any(!is.finite(lk$score))) break
      score <- lk$score + drop(P %*% beta)
      K <- lk$K + P
      R <- tryCatch(chol(K), error = function(e) NULL)
      if(is.null(R)) {
        shift <- max(1e-06, abs(min(eigen(K, symmetric = TRUE,
          only.values = TRUE)$values)) * 1.1)
        R <- tryCatch(chol(K + diag(shift, p)), error = function(e) NULL)
        if(is.null(R)) break
      }
      step <- drop(backsolve(R, forwardsolve(t(R), score)))
      if(max(abs(step)) < 1e-09 * max(1, abs(beta))) {
        good <- TRUE
        break
      }
      f0 <- lk$value + 0.5 * drop(crossprod(beta, P %*% beta))
      alpha <- 1
      for(i in seq_len(25L)) {
        trial <- beta - alpha * step
        f1 <- likelihood(trial, FALSE) +
          0.5 * drop(crossprod(trial, P %*% trial))
        if(is.finite(f1) && f1 <= f0 + 1e-10) break
        alpha <- alpha / 2
      }
      if(!is.finite(f1) || f1 > f0 + 1e-10) break
      beta <- trial
      if(max(abs(alpha * step)) < 1e-09 * max(1, abs(beta))) {
        good <- TRUE
        break
      }
    }
    value <- Inf
    V <- NULL
    if(good) {
      lk <- likelihood(beta)
      K <- lk$K + P
      R <- if(all(is.finite(K)))
        tryCatch(chol(K), error = function(e) NULL) else NULL
      if(!is.null(R)) {
        V <- chol2inv(R)
        value <- unname(lk$value + 0.5 *
          (drop(crossprod(beta, P %*% beta)) +
            2 * sum(log(diag(R))) - logdetP))
      }
    }
    if(is.finite(value)) state$beta <- beta
    state$last.rho <- rho
    state$last.value <- value
    state$last.beta <- if(is.finite(value)) beta else NULL
    state$last.V <- V
    value
  }

  if(!is.finite(profile(rho0))) {
    reason <- "initial joint REML evaluation failed"
    ## dBCT switches to its BCCG limit above tau = 1e6. An RS fit
    ## approaching that limit has no usable curvature in the tau direction.
    if(identical(kind, "BCT") && "tau" %in% np &&
        identical(unname(family$links["tau"]), "log") &&
        all(is.finite(eta0$tau)) &&
        diff(range(eta0$tau)) < 1e-10) {
      tau <- exp(eta0$tau[1L])
      if(is.finite(tau) && tau >= 1e4 && tau < 5e5) {
        trial <- eta0
        trial$tau <- trial$tau + log(2)
        next.loglik <- tryCatch(sum(w * family$pdf(
          par = family$map2par(trial), y = y, log = TRUE)),
          error = function(e) NA_real_)
        if(is.finite(next.loglik) && next.loglik > -check + 1e-08)
          reason <- sprintf(paste0("BCT tau approaches its large-tau limit ",
            "(RS tau = %.0f; log likelihood rises when tau is doubled); ",
            "a regular finite joint mode was not found"), tau)
      }
    }
    return(unavailable(reason))
  }
  finite_gradient <- function(rho) {
    h <- 0.01
    g <- numeric(length(rho))
    f0 <- profile(rho)
    for(i in seq_along(rho)) {
      r <- rho
      r[i] <- min(upper[i], rho[i] + h)
      fp <- profile(r)
      r[i] <- max(lower[i], rho[i] - h)
      fm <- profile(r)
      if(is.finite(fp) && is.finite(fm)) {
        g[i] <- (fp - fm) /
          (min(upper[i], rho[i] + h) - max(lower[i], rho[i] - h))
      } else if(is.finite(fp) && is.finite(f0)) {
        g[i] <- (fp - f0) / (min(upper[i], rho[i] + h) - rho[i])
      } else if(is.finite(fm) && is.finite(f0)) {
        g[i] <- (f0 - fm) / (rho[i] -
          max(lower[i], rho[i] - h))
      } else {
        g[i] <- 0
      }
    }
    g
  }
  gradient.mode <- control$jr.gradient
  if(is.null(gradient.mode))
    gradient.mode <- if(length(rho0) < 3L) "finite" else "implicit"
  if(!gradient.mode %in% c("finite", "implicit", "hybrid"))
    stop("jr.gradient must be 'finite', 'implicit', or 'hybrid'.")
  gradient.fallbacks <- 0L
  implicit_gradient <- function(rho) {
    profile(rho)
    beta <- state$last.beta
    V <- state$last.V
    if(is.null(beta) || is.null(V))
      stop("direct coefficient mode is unavailable")
    g <- numeric(length(rho))
    for(z in blocks) {
      if(!length(z$ids)) next
      lambda <- exp(rho[z$ids])
      Pr <- NULL
      if(!is.null(z$basis) && length(z$S) > 1L) {
        Pr <- matrix(0, z$basis$rank, z$basis$rank)
        for(i in seq_along(z$S))
          Pr <- Pr + lambda[i] * z$basis$penalties[[i]]
        Pr <- chol2inv(chol(Pr))
      }
      for(i in seq_along(z$ids)) {
        Pi <- lambda[i] * z$S[[i]]
        u <- numeric(p)
        u[z$index] <- drop(Pi %*% beta[z$index])
        db <- -drop(V %*% u)
        deta <- lapply(np, function(j)
          drop(X[[j]] %*% db[index[[j]]]))
        scale <- max(abs(unlist(deta, use.names = FALSE)))
        h <- if(scale > 0) min(0.05, 0.01 / scale) else 0.05
        Hp <- likelihood(beta + h * db)$K
        Hm <- likelihood(beta - h * db)$K
        if(any(!is.finite(Hp)) || any(!is.finite(Hm)))
          stop("direct likelihood curvature is unavailable")
        penalty.trace <- if(is.null(z$basis)) 0 else
          if(length(z$S) == 1L) z$basis$rank else
            sum(Pr * (lambda[i] * z$basis$penalties[[i]]))
        g[z$ids[i]] <- 0.5 * (sum(beta[z$index] * u[z$index]) +
          sum(V[z$index, z$index, drop = FALSE] * Pi) -
          penalty.trace + sum(V * (Hp - Hm)) / (2 * h))
      }
    }
    if(any(!is.finite(g))) stop("joint REML gradient is not finite")
    g
  }
  outer_gradient <- if(gradient.mode == "finite") finite_gradient else
    function(rho) {
      g <- tryCatch(implicit_gradient(rho), error = function(e) NULL)
      if(!is.null(g)) return(g)
      gradient.fallbacks <<- gradient.fallbacks + 1L
      finite_gradient(rho)
    }
  opt <- tryCatch(nlminb(rho0, profile, gradient = outer_gradient,
    lower = lower, upper = upper,
    control = list(iter.max = if(is.null(control$jr.outer.maxit)) 30L else
      control$jr.outer.maxit, eval.max = if(is.null(control$jr.outer.maxeval))
      150L else control$jr.outer.maxeval, rel.tol = 1e-08)),
    error = function(e) e)
  if(inherits(opt, "error"))
    return(unavailable(conditionMessage(opt)))
  rho <- opt$par
  value <- profile(rho)
  if(!is.finite(value))
    return(unavailable("optimized joint REML evaluation failed"))
  beta <- state$last.beta
  Vmode <- state$last.V
  dimnames(Vmode) <- list(names(co), names(co))
  Pmode <- matrix(0, p, p)
  for(z in blocks) {
    lambda <- if(length(z$ids)) exp(rho[z$ids]) else z$lambda
    for(i in seq_along(z$S))
      Pmode[z$index, z$index] <- Pmode[z$index, z$index] +
        lambda[i] * z$S[[i]]
  }
  lk <- likelihood(beta)
  mode.residual <- max(abs(drop(Vmode %*%
    (lk$score + drop(Pmode %*% beta))))) / max(1, abs(beta))
  residual <- mode.residual
  if(!is.finite(mode.residual) || mode.residual > 1e-04)
    return(unavailable("direct coefficient mode could not be verified"))

  ## Update the stored fit from the direct joint mode. An RS refit would
  ## solve different family working equations for some distributions.
  fit <- fit0
  fit$family <- family
  for(j in np) {
    linear <- fit$coefficients[[j]]
    if(length(linear)) {
      ii <- match(paste0(j, ".p.", names(linear)), names(co))
      if(anyNA(ii))
        stop("JR could not match the stored linear coefficients.")
      fit$coefficients[[j]][] <- beta[ii]
      fit$fitted.linear[[j]]$coefficients[] <- beta[ii]
      fit$fitted.linear[[j]]$fitted.values <-
        drop(X[[j]][, match(ii, index[[j]]), drop = FALSE] %*% beta[ii])
      fit$fitted.linear[[j]]$vcov <- Vmode[ii, ii, drop = FALSE]
    }
  }
  for(z in blocks) {
    j <- z$parameter
    k <- z$term
    ii <- z$index
    lambda <- if(length(z$ids)) exp(rho[z$ids]) else z$lambda
    term <- fit$fitted.specials[[j]][[k]]
    term$coefficients[] <- beta[ii]
    term$lambdas <- lambda
    term$fitted.values <-
      drop(X[[j]][, match(ii, index[[j]]), drop = FALSE] %*% beta[ii])
    term$vcov <- Vmode[ii, ii, drop = FALSE]
    edf <- length(ii) - sum(term$vcov * Pmode[ii, ii, drop = FALSE])
    term$edf <- edf
    term$df <- n - edf
    fit$fitted.specials[[j]][[k]] <- term
  }
  eta <- lapply(np, function(j)
    drop(X[[j]] %*% beta[index[[j]]]) + off[[j]])
  names(eta) <- np
  fit$fitted.values <- as.data.frame(eta)
  fit$logLik <- -likelihood(beta, FALSE)
  fit$deviance <- -2 * fit$logLik
  fit$dev.reduction <- (fit$null.deviance - fit$deviance) /
    fit$null.deviance
  fit$df <- get_df(fit)

  ## Profile curvature and covariance-factor derivatives use the same solver.
  m <- length(rho)
  h <- rep(0.01, m)
  gradient <- rep(NA_real_, m)
  H <- matrix(NA_real_, m, m)
  dC <- vector("list", m)
  covariance_factor <- function(V) {
    R <- chol(chol2inv(chol(V)))
    backsolve(R, diag(nrow(V)))
  }
  if(all(rho - h > lower & rho + h < upper)) {
    for(i in seq_len(m)) {
      r <- rho
      r[i] <- rho[i] + h[i]
      fp <- profile(r)
      Vp <- state$last.V
      r[i] <- rho[i] - h[i]
      fm <- profile(r)
      Vm <- state$last.V
      gradient[i] <- (fp - fm) / (2 * h[i])
      H[i, i] <- (fp - 2 * value + fm) / h[i]^2
      if(!is.null(Vp) && !is.null(Vm)) {
        Cp <- tryCatch(covariance_factor(Vp), error = function(e) NULL)
        Cm <- tryCatch(covariance_factor(Vm), error = function(e) NULL)
        if(!is.null(Cp) && !is.null(Cm))
          dC[[i]] <- (Cp - Cm) / (2 * h[i])
      }
    }
    if(m > 1L) for(i in seq_len(m - 1L)) for(j in (i + 1L):m) {
      v <- numeric(4L)
      signs <- rbind(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
      for(k in seq_len(4L)) {
        r <- rho
        r[i] <- r[i] + signs[k, 1L] * h[i]
        r[j] <- r[j] + signs[k, 2L] * h[j]
        v[k] <- profile(r)
      }
      H[i, j] <- H[j, i] <-
        (v[1L] - v[2L] - v[3L] + v[4L]) / (4 * h[i] * h[j])
    }
  }
  covariance <- factor.root <- NULL
  curvature.rank <- 0L
  curvature.cutoff <- NA_real_
  curvature.values <- rep(NA_real_, m)
  stationary <- FALSE
  if(all(is.finite(H)) && all(is.finite(gradient))) {
    ev <- eigen(0.5 * (H + t(H)), symmetric = TRUE)
    curvature.values <- ev$values
    curvature.cutoff <- max(1e-04 * max(abs(ev$values)),
      100 * .Machine$double.eps * max(1, abs(value)) / min(h)^2)
    positive <- ev$values > curvature.cutoff
    curvature.rank <- sum(positive)
    if(any(positive) && !any(ev$values < -curvature.cutoff)) {
      Q <- ev$vectors[, positive, drop = FALSE]
      Q <- sweep(Q, 2L, sqrt(ev$values[positive]), "/")
      covariance <- tcrossprod(Q)
      step <- drop(covariance %*% gradient)
      flat.gradient <- gradient -
        drop(ev$vectors[, positive, drop = FALSE] %*%
          crossprod(ev$vectors[, positive, drop = FALSE], gradient))
      stationary <- max(abs(step)) < 0.05 &&
        sum(gradient * step) < 0.04 &&
        max(abs(flat.gradient)) < 1e-05 &&
        all(vapply(dC, function(z) !is.null(z) && all(is.finite(z)),
          logical(1L)))
      factor.variance <- rep(if(length(np) > 1L) 50 else 10, m)
      factor.variance[positive] <- 1 / ev$values[positive]
      factor.root <- t(sweep(ev$vectors, 2L,
        sqrt(factor.variance), "*"))
    }
  }
  if(!stationary) {
    covariance <- NULL
    factor.root <- NULL
    dC <- NULL
  }
  names(rho) <- names(rho0)
  names(gradient) <- names(rho0)
  if(!is.null(dC)) names(dC) <- names(rho0)
  dimnames(H) <- list(names(rho0), names(rho0))
  if(!is.null(covariance)) dimnames(covariance) <- dimnames(H)
  if(!is.null(factor.root)) colnames(factor.root) <- names(rho0)
  fit$jr <- list(initial = initial, rho = rho, criterion = value,
    gradient = gradient, hessian = H, covariance = covariance,
    coefficient.covariance = Vmode,
    mode.coefficients = setNames(beta, names(co)),
    factor.root = factor.root, curvature.rank = curvature.rank,
    curvature.cutoff = curvature.cutoff, curvature.values = curvature.values,
    factor.derivative = dC, fixed.point.residual = residual,
    joint.mode.residual = mode.residual,
    score.source = if(analytic) setNames(rep("closed-form", length(np)), np)
      else ifelse(numerical.score, "density", "family"),
    converged = opt$convergence == 0L && stationary &&
      is.finite(residual) && residual < 1e-05 &&
      is.finite(mode.residual) && mode.residual < 1e-04,
    evaluations = state$evaluations, gradient.mode = gradient.mode,
    gradient.fallbacks = gradient.fallbacks, outer = opt)
  fit$control <- control
  fit$control$.jr.initial.fit <- NULL
  fit$converged <- isTRUE(fit$converged) && fit$jr$converged
  if(!fit$jr$converged)
    return(unavailable(
      "joint REML curvature or coefficient fixed point could not be verified"))
  if(isTRUE(control$trace))
    cat("GAMLSS-JR: joint REML = ", signif(value, 8L),
      " evaluations = ", state$evaluations,
      " converged = ", fit$jr$converged, "\n", sep = "")
  fit
}


## Continue fast joint REML estimation from a complete fitted model.
jr <- function(object, trace = TRUE, ...)
{
  if(!inherits(object, "gamlss2") || inherits(object, "bamlss2"))
    stop("jr requires a fitted gamlss2 model.")
  if(is.null(object$x) || is.null(object$y) ||
      is.null(object$specials) || is.null(object$xterms) ||
      is.null(object$sterms))
    stop("jr requires a fit retaining x, y, and the smooth setup.")
  if(any(vapply(object$specials, function(z)
      inherits(z, "mgcv.smooth") && is.null(z$X), logical(1L))))
    stop("jr requires the stored smooth design matrices.")

  control <- object$control
  control$optimizer <- JR
  control$trace <- trace
  extra <- list(...)
  if(length(extra)) {
    if(is.null(names(extra)) || any(!nzchar(names(extra))))
      stop("jr control arguments must be named.")
    for(k in names(extra)) control[[k]] <- extra[[k]]
  }
  control$.jr.initial.fit <- object
  tstart <- proc.time()
  fit <- JR(x = object$x, y = object$y, specials = object$specials,
    family = object$family, offsets = object$offsets,
    weights = object$weights,
    start = coef(object, full = TRUE, lambdas = TRUE, dropall = FALSE),
    xterms = object$xterms, sterms = object$sterms, control = control)
  elapsed <- as.numeric((proc.time() - tstart)["elapsed"])
  if(is.null(fit[["jr", exact = TRUE]])) {
    object$reml.failure <- fit$reml.failure
    object$reml.failure$elapsed <- elapsed
    return(object)
  }

  ## Keep the original model call and setup for prediction and plotting.
  object$reml.failure <- NULL
  for(k in names(fit)) object[[k]] <- fit[[k]]
  object$df <- get_df(object)
  object$elapsed <- object$elapsed + elapsed
  object$jr$elapsed <- elapsed
  object$jr.call <- match.call()
  if(!is.null(object$call))
    object$call$optimizer <- quote(JR)
  if(!is.null(object$model))
    object$results <- results(object, data = object$model, interval = "local")
  object
}
