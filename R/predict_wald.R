wald_information <- function(object)
{
  if(inherits(object, "bamlss2"))
    stop("Wald intervals require an ML fit; use FUN for posterior intervals.")
  if(isTRUE(object$control$ridge))
    stop("Wald intervals do not yet support automatically selected linear ridge penalties.")
  if(is.null(object$y) || is.null(object$fitted.values))
    stop("Wald intervals require a fit retaining y and fitted values (light = FALSE, y = TRUE).")

  family <- object$family
  parameters <- family$names
  eta <- as.list(object$fitted.values[parameters])
  n <- length(eta[[1L]])
  y <- object$y
  X <- model.matrix(object)
  blocks <- designs <- indices <- setNames(vector("list", length(parameters)), parameters)
  penalties <- list()
  dimension <- 0L
  for(j in parameters) {
    fixed <- object$control$fixed[[j]]
    fixed <- length(fixed) == 1L && !is.na(fixed) && as.logical(fixed)
    add_block <- function(B, beta, label = NULL, S = NULL, lambda = NULL) {
      if(anyNA(beta) || any(!is.finite(beta)))
        stop("Wald intervals require estimable, finite coefficients; parameter '", j, "'.")
      if(!is.matrix(B) || nrow(B) != n || ncol(B) != length(beta) ||
          any(!is.finite(B)))
        stop("Invalid prediction design for parameter '", j, "'.")
      ind <- if(fixed) integer() else dimension + seq_along(beta)
      blocks[[j]][[length(blocks[[j]]) + 1L]] <<- list(
        label = label, names = names(beta), index = ind, size = length(beta))
      if(!fixed) {
        dimension <<- dimension + length(beta)
        designs[[j]] <<- cbind(designs[[j]], B)
        indices[[j]] <<- c(indices[[j]], ind)
        if(length(S)) {
          if(!is.list(S)) S <- list(S)
          if(is.null(lambda) || !length(lambda))
            stop("Missing fitted smoothing parameters for term '", label, "'.")
          lambda <- rep(lambda, length.out = length(S))
          P <- matrix(0, length(beta), length(beta))
          for(k in seq_along(S)) {
            if(!identical(dim(S[[k]]), dim(P)) || !is.finite(lambda[k]))
              stop("Invalid smoothing penalty for term '", label, "'.")
            P <- P + lambda[k] * S[[k]]
          }
          penalties[[length(penalties) + 1L]] <<- list(index = ind, P = P)
        }
      }
    }
    beta <- object$coefficients[[j]]
    if(length(beta)) add_block(X[, names(beta), drop = FALSE], beta)
    for(term in names(object$fitted.specials[[j]])) {
      fitted <- object$fitted.specials[[j]][[term]]
      smooth <- object$specials[[term]]
      if(!inherits(smooth, "mgcv.smooth") || !is.null(smooth$special.wfit))
        stop("Wald intervals currently support linear terms and standard mgcv smooths; unsupported term '", term, "'.")
      if(isTRUE(object$control$termselect) || isTRUE(smooth$control$termselect))
        stop("Wald intervals do not yet support adaptive term-selection penalties.")
      B <- smooth$X
      if(!is.null(smooth$binning))
        B <- B[smooth$binning$match.index, , drop = FALSE]
      add_block(B, fitted$coefficients, term, smooth[["S", exact = TRUE]], fitted$lambdas)
    }
  }

  H <- matrix(0, dimension, dimension)
  weights <- object$weights
  if(is.null(weights)) weights <- rep(1, n)
  if(length(weights) != n || any(!is.finite(weights)) || any(weights < 0))
    stop("Wald intervals require finite, nonnegative observation weights.")
  active <- parameters[lengths(indices) > 0L]
  step <- .Machine$double.eps^(1/3)
  for(b in active) {
    h <- step * pmax(1, abs(eta[[b]]))
    upper <- lower <- eta
    upper[[b]] <- eta[[b]] + h
    lower[[b]] <- eta[[b]] - h
    upper <- family$map2par(upper)
    lower <- family$map2par(lower)
    for(a in active) {
      score <- family$score[[a]]
      if(!is.function(score)) stop("Missing linked score for parameter '", a, "'.")
      su <- as.numeric(score(par = upper, y = y))
      sl <- as.numeric(score(par = lower, y = y))
      if(length(su) != n || length(sl) != n)
        stop("Wald intervals require one score per observation.")
      curvature <- -(su - sl) / (2 * h)
      curvature[weights == 0] <- 0
      curvature <- curvature * weights
      if(any(!is.finite(curvature)))
        stop("Non-finite score curvature; analytical Wald intervals are unavailable.")
      H[indices[[a]], indices[[b]]] <- crossprod(designs[[a]], designs[[b]] * curvature)
    }
  }
  H <- (H + t(H)) / 2
  for(penalty in penalties)
    H[penalty$index, penalty$index] <- H[penalty$index, penalty$index] + penalty$P
  R <- if(dimension) tryCatch(chol(H), error = function(e) NULL) else H
  if(is.null(R))
    stop("The joint information matrix is not positive definite; Wald intervals are unavailable. Check convergence/identifiability or use the existing simulation method.")
  structure(list(object = object, R = R, blocks = blocks, indices = indices),
    class = "gamlss2.interval.cache")
}

## Only solve for the requested prediction variances, in bounded row chunks.
wald_variance <- function(A, R)
{
  n <- nrow(A)
  variance <- numeric(n)
  if(!n || !ncol(A)) return(variance)
  for(first in seq.int(1L, n, by = 1024L)) {
    rows <- seq.int(first, min(n, first + 1023L))
    B <- backsolve(R, t(A[rows, , drop = FALSE]), transpose = TRUE)
    variance[rows] <- colSums(B * B)
  }
  variance
}

print.gamlss2.interval.cache <- function(x, ...)
{
  cat("Reusable Wald information factor:", ncol(x$R), "coefficients\n")
  invisible(x)
}

predict_wald <- function(object, model, newdata, type, terms, drop, dots,
  level, cache)
{
  if(!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1)
    stop("'level' must be a number strictly between zero and one.")
  if(any(c("FUN", "R", "seed", "burnin") %in% names(dots)))
    stop("Do not combine interval = 'wald' with FUN, R, seed, or burnin; these select simulation summaries.")
  if(inherits(object, "bamlss2"))
    stop("Wald intervals require an ML fit; use FUN for posterior intervals.")
  if(is.null(cache)) {
    cache <- wald_information(object)
  } else if(!inherits(cache, "gamlss2.interval.cache") ||
      !identical(cache$object, object)) {
    stop("'interval.cache' must come from the same, unchanged fitted model.")
  }
  family <- object$family
  parameters <- family$names
  if(is.null(model)) model <- dots$what
  if(is.null(model)) model <- dots$parameter
  if(is.null(model)) model <- parameters
  if(!is.character(model)) model <- parameters[model]
  model <- parameters[pmatch(model, parameters)]
  if(!length(model) || anyNA(model) || anyDuplicated(model))
    stop("Unknown or duplicated prediction parameter.")
  if(type == "response" && !setequal(model, parameters))
    stop("Response intervals require all distributional parameters.")

  ## Ordinary predictions remain the source of the point estimates and their
  ## output shape. This call never requests draws or standard errors.
  if(is.null(dots$no_weights)) dots$no_weights <- TRUE
  point <- function(type, model, drop) do.call(predict.gamlss2,
    c(list(object = object, model = model, newdata = newdata, type = type,
      terms = terms, drop = drop), dots))
  fit <- point(type, model, drop)
  mf.args <- c(list(formula = object), dots)
  if(!is.null(newdata)) mf.args$data <- newdata
  mf <- do.call(model.frame.gamlss2, mf.args)
  X <- model.matrix(object, data = mf)
  n <- nrow(X)
  p <- ncol(cache$R)
  matrices <- term.matrices <- setNames(vector("list", length(parameters)), parameters)
  terms <- gsub(" ", "", terms)
  if(!length(terms)) terms <- NULL
  for(j in parameters) {
    labels <- c(object$xterms[[j]], object$sterms[[j]])
    selected <- if(is.null(terms)) labels else if(isTRUE(dots$nogrep)) terms else
      grep2(terms, labels, fixed = TRUE, value = TRUE)
    selected <- unique(selected)
    wanted <- function(labels) {
      if(isTRUE(dots$nogrep)) labels[labels %in% selected] else
        unique(unlist(lapply(selected, grep, x = labels, fixed = TRUE, value = TRUE)))
    }
    linear <- wanted(object$xterms[[j]])
    smooths <- wanted(object$sterms[[j]])
    A <- matrix(0, n, p)
    term.matrices[[j]] <- list()
    for(block in cache$blocks[[j]]) {
      if(is.null(block$label)) {
        use <- which(block$names %in% linear)
        if(length(block$index) && length(use))
          A[, block$index[use]] <- X[, block$names[use], drop = FALSE]
        for(k in use) {
          term.matrices[[j]][[block$names[k]]] <- list(
            index = if(length(block$index)) block$index[k] else integer(),
            X = X[, block$names[k], drop = FALSE])
        }
      } else if(block$label %in% smooths) {
        B <- PredictMat(object$specials[[block$label]], mf, n = n)
        if(length(block$index))
          A[, block$index] <- B
        term.matrices[[j]][[block$label]] <- list(index = block$index, X = B)
      }
    }
    matrices[[j]] <- A
  }
  critical <- qnorm((1 + level) / 2)
  se <- lower <- upper <- fit
  if(type == "terms") {
    for(j in model) {
      values <- if(length(model) == 1L && drop) fit else fit[[j]]
      sj <- lj <- uj <- values
      for(term in colnames(values)) {
        block <- term.matrices[[j]][[term]]
        A <- matrix(0, n, p)
        if(length(block$index)) A[, block$index] <- block$X
        sj[, term] <- sqrt(wald_variance(A, cache$R))
        lj[, term] <- values[, term] - critical * sj[, term]
        uj[, term] <- values[, term] + critical * sj[, term]
      }
      if(length(model) == 1L && drop) {
        se <- sj; lower <- lj; upper <- uj
      } else {
        se[[j]] <- sj; lower[[j]] <- lj; upper[[j]] <- uj
      }
    }
  } else {
    eta <- if(type != "link") as.list(point("link", parameters, FALSE)) else NULL
    fm <- family$mean
    if(is.null(fm)) fm <- if(!is.null(family$q))
      function(par) family$quantile(par = par, 0.5) else function(par) par[[1L]]
    targets <- if(type == "response") "response" else model
    ## Small predictor-space Jacobian for transformed parameters or means.
    derivatives <- setNames(vector("list", length(targets)), targets)
    if(type != "link") {
      for(j in parameters) {
        h <- .Machine$double.eps^(1/3) * pmax(1, abs(eta[[j]]))
        hi <- lo <- eta
        hi[[j]] <- eta[[j]] + h; lo[[j]] <- eta[[j]] - h
        hi <- family$map2par(hi); lo <- family$map2par(lo)
        for(target in targets)
          derivatives[[target]][[j]] <- if(type == "response")
            (fm(hi) - fm(lo)) / (2 * h) else (hi[[target]] - lo[[target]]) / (2 * h)
      }
    }
    for(target in targets) {
      values <- if(type == "response" || (length(model) == 1L && drop)) fit else fit[[target]]
      A <- if(type == "link") matrices[[target]] else matrix(0, n, p)
      if(type != "link") for(j in parameters)
        A <- A + matrices[[j]] * as.numeric(derivatives[[target]][[j]])
      sj <- values
      sj[] <- sqrt(wald_variance(A, cache$R))
      lj <- values - critical * sj; uj <- values + critical * sj
      if(type == "parameter") {
        link <- family$links[[target]]
        others <- setdiff(parameters, target)
        separable <- all(vapply(others, function(j)
          all(derivatives[[target]][[j]] == 0), logical(1L)))
        if(is.character(link) && length(link) == 1L &&
            link %in% c("identity", "log", "logit", "probit", "cloglog", "cauchit") && separable) {
          link <- make.link2(link)
          if(isTRUE(all.equal(as.numeric(link$linkinv(eta[[target]])),
              as.numeric(values), tolerance = 1e-12))) {
            s <- sqrt(wald_variance(matrices[[target]], cache$R))
            lj <- link$linkinv(eta[[target]] - critical * s)
            uj <- link$linkinv(eta[[target]] + critical * s)
          }
        }
      }
      if(type == "response" || (length(model) == 1L && drop)) {
        se <- sj; lower <- lj; upper <- uj
      } else {
        se[[target]] <- sj; lower[[target]] <- lj; upper[[target]] <- uj
      }
    }
  }
  structure(list(fit = fit, se.fit = se, lower = lower, upper = upper),
    level = level, interval = "wald", interval.cache = cache)
}
