## auxiliary function for checking whether dedicated S3 methods are available
hasS3method <- function (method, classes) {
  any(sapply(classes, function(cls) {
    tryCatch(is.function(getS3method(method, class = cls)), error = function(e) FALSE)
  }))
}

## S3 method for generating a gamlss2.family from a distribution object (from distributions3)
family.distribution <- function(object, links, score = TRUE, hessian = FALSE, update = FALSE, ...) {
  ## distributions3 class
  d <- setdiff(class(object), "distribution")
  create_distribution <- function(par) structure(par, class = class(object))

  ## score and hessian methods available?
  has_score <- hasS3method("score", d)
  has_hessian <- hasS3method("hessian", d)
  
  distributions3_family(d[1L],
    links = links,
    score = score && has_score,
    hessian = hessian && has_hessian,
    update = update && has_score && has_hessian,
    create_distribution = create_distribution,
    ...)
}

## S3 methods for generating a gamlss2.family from gamlss.dist objects
## (either a gamlss.family or a GAMLSS distribution)
family.GAMLSS <- function(object, ...) {
  object <- gamlss.dist:::.gamlss_family(object)
  NextMethod()
}

family.gamlss.family <- function(object, score = TRUE, hessian = FALSE, update = FALSE, ...) {
  ## links specifications
  links <- unlist(object[endsWith(names(object), "link")])
  names(links) <- gsub(".link", "", names(links), fixed = TRUE)
  
  ## initialize functions
  initialize <- lapply(names(links), function(p) {
    function(y, ...) as.numeric(eval(object[[paste0(p, ".initial")]]))
  })
  names(initialize) <- names(links)
  
  ## distributions3 class
  create_distribution <- function(par) structure(par, class = c("GAMLSS", "distribution"), family = object$family)

  distributions3_family("GAMLSS", links = links, score = score, hessian = hessian, update = update,
    create_distribution = create_distribution, type = tolower(object$type), initialize = initialize, ...)
}

## workhorse function
distributions3_family <- function(distribution, links, score = TRUE, hessian = FALSE, update = FALSE, create_distribution = "structure",
  type = NULL, valid.response = NULL, initialize = NULL, ...) {

  ## distribution must be a character,
  ## corresponding to a distribution generator
  stopifnot(is.character(distribution))
  args <- argsAnywhere(distribution)
  stopifnot(is.function(args))
  nams <- names(formals(args))

  ## default strategies for creating distributions3 objects
  if (is.character(create_distribution)) create_distribution <- match.arg(create_distribution, c("structure", "do.call"))

  ## check for extra arguments
  xtra <- list(...)
  if (length(xtra) > 0L) {
    nams <- setdiff(nams, names(xtra))
    if (is.character(create_distribution) && create_distribution == "structure") {
      warning("'create_distribution' must be a function or set to 'do.call' for distributions with extra arguments")
      create_distribution <- "do.call"
    }
  }

  ## links must be provided as named vector or list (with matching names)
  links <- as.list(links)
  stopifnot(is.list(links), !is.null(names(links)))
  nams <- intersect(nams, names(links))
  if (length(nams) < 1L) warning("'distribution' and 'links' do not seem to have any parameter names in common")
  links <- links[nams]
  for (n in nams) {
    if (is.character(links[[n]])) links[[n]] <- make.link2(links[[n]])
    stopifnot(inherits(links[[n]], "link-glm"))
  }

  ## create distributions3 object from parameter list
  d3 <- if (is.function(create_distribution)) {
    create_distribution
  } else if (is.character(create_distribution) && create_distribution == "structure") {
      function(par) structure(par, class = c(distribution, "distribution"))
  } else if (is.character(create_distribution) && create_distribution == "do.call") {
      function(par) do.call(distribution, c(par, xtra))
  } else {
    stop("'create_distribution' must either be a function(par) or a character string")
  }
  
  ## set up family list
  rval <- list(
    "family"   = distribution,
    "names"    = nams,
    "links"    = links,
    "log_likelihood"   = function(par, y, ...) sum(distributions3::log_pdf(d3(par), y, elementwise = TRUE, ...)),
    "mu"       = function(par, ...) mean(d3(par), ...),
    "pdf"      = function(par, y, log = FALSE) distributions3::pdf(d3(par), y, elementwise = TRUE, log = log),
    "cdf"      = function(par, y, ...) distributions3::cdf(d3(par), y, elementwise = TRUE, ...),
    "random"   = function(par, n) distributions3::random(d3(par), n),
    "quantile" = function(par, p) quantile(d3(par), p, elementwise = TRUE),
    "crps"     = function(par, y, ...) sum(scoringRules::crps(d3(par), y, elementwise = TRUE, ...)),
    "mean"     = function(par) mean(d3(par)),
    "variance" = function(par) distributions3::variance(d3(par)),
    "skewness" = function(par) distributions3::skewness(d3(par)),
    "kurtosis" = function(par) distributions3::kurtosis(d3(par)),
    "create_distribution" = d3
  )

  ## add score function if desired and available
  has_score <- hasS3method("score", distribution)
  if (score && has_score) rval$score <- structure(lapply(nams, function(n) { ## FIXME: score(par, y, which, ...)
    function(par, y, ...) {
      lnk <- links[[n]]
      eta <- lnk$linkfun(par[[n]])
      score(d3(par), y, which = n, ...) * lnk$mu.eta(eta)
    }
  }), names = nams)

  ## add hessian function if desired and available
  has_hessian <- hasS3method("hessian", distribution)
  if (hessian && has_hessian) rval$hess <- structure(lapply(nams, function(n) { ## FIXME: hessian(par, y, which, ...)
    function(par, y, ...) {
      lnk <- links[[n]]
      eta <- lnk$linkfun(par[[n]])
      par <- d3(par)
      s <- score(par, y, which = n, ...)
      h <- hessian(par, y, which = n, ...)
      h * lnk$mu.eta(eta)^2 + s * lnk$mu.eta2(eta)
    }
  }), names = nams)

  ## add Newton/Fisher update function if desired and available
  if (update && has_score && has_hessian) rval$z_weights <- function(y, eta, peta, j) { ## FIXME: update(par, y, eta, which) -> list(eta, weights)
    lnk <- links[[j]]
    mu.eta <- lnk$mu.eta(eta)
    par <- d3(peta)

    s <- score(par, y, which = j)
    score <- s * mu.eta
    score <- gamlss2:::deriv_checks(score, is.weight = FALSE)

    h <- hessian(par, y, which = j)
    hess <- h * mu.eta^2 + s * lnk$mu.eta2(eta)
    hess <- gamlss2:::deriv_checks(hess, is.weight = TRUE)

    list(
      z = eta + 1/hess * score,
      weights = hess
    )
  }

  ## family elements that are difficult to infer automatically
  rval$valid.response <- if (is.null(valid.response)) function(x) TRUE else valid.response
  if (!is.null(type)) rval$type <- match.arg(type, c("continuous", "discrete")) ## FIXME: further types allowed?
  if (!is.null(initialize) && is.list(initialize) && all(names(initialize) %in% nams)) rval$initialize <- initialize

  ## return classed list
  class(rval) <- "gamlss2.family"
  return(rval)
}



## -----------------------------------------------------------------------------

## the following code is intended for distributions3 but has to wait for mergin a major PR there

## score/hessian generics
score <- function(d, ...) UseMethod("score")
hessian <- function(d, ...) UseMethod("hessian")

## fallback methods based on numeric differentiation
score.distribution <- function(d, x, which = NULL, drop = TRUE, eps = .Machine$double.eps^(1/3), ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## available and selected parameters
  p <- names(unclass(d))
  if (is.null(which)) which <- p
  which <- match.arg(which, p, several.ok = TRUE)

  ## compute scores
  scr <- function(par) {
    d1 <- d2 <- d
    d1[[par]] <- d1[[par]] + eps
    d2[[par]] <- d2[[par]] - eps
    (log_pdf(d1, x) - log_pdf(d2, x)) / (2 * eps)
  }

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    s <- setNames(scr(which), names(d))
  } else {
    s <- lapply(which, scr)
    s <- do.call("cbind", s)
    dimnames(s) <- list(names(d), which)
  }
  return(s)
}

hessian.distribution <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, eps = .Machine$double.eps^(1/4), ...) {
  ## numeric differentiation yields observed hessian only
  if (!identical(expected, FALSE)) stop("only the observed hessian is available")
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## available and selected parameters
  p <- names(unclass(d))
  pp <- outer(p, p, paste, sep = ":")
  diag(pp) <- p
  p <- setNames(
    c(diag(pp), pp[upper.tri(pp)], pp[upper.tri(pp)]),
    c(diag(pp), pp[upper.tri(pp)], pp[lower.tri(pp)])
  )[pp]
  if (is.null(which)) which <- names(p)
  
  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ## compute scores
  hess <- function(par) {
    par <- strsplit(par, ":", fixed = TRUE)[[1L]]
    par <- rep_len(par, 2L)
    d1 <- d2 <- d
    d1[[par[2L]]] <- d1[[par[2L]]] + eps
    d2[[par[2L]]] <- d2[[par[2L]]] - eps
    (score(d1, x, which = par[1L]) - score(d2, x, which = par[1L])) / (2 * eps)
  }

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    h <- setNames(hess(w), names(d))
  } else {
    h <- lapply(w, hess)
    h <- do.call("cbind", h)
    dimnames(h) <- list(names(d), w)
    if (!identical(w, which)) h <- h[, p[which], drop = FALSE]
    colnames(h) <- which
  }
  return(h)
}


## Normal methods for score/hessian
score.Normal <- function(d, x, which = NULL, drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## available and selected parameters
  p <- c("mu", "sigma")
  if (is.null(which)) which <- p
  which <- match.arg(which, p, several.ok = TRUE)

  ## compute scores
  scr <- function(par) switch(par,
    "mu"    = (x - d$mu)/(d$sigma^2),
    "sigma" = (x - d$mu)^2/(d$sigma^3) - 1/d$sigma)

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    s <- setNames(scr(which), names(d))
  } else {
    s <- lapply(which, scr)
    s <- do.call("cbind", s)
    dimnames(s) <- list(names(d), which)
  }
  return(s)
}

hessian.Normal <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## available and selected parameters/combinations and mappings for symmetries
  p <- c("mu" = "mu", "sigma:mu" = "mu:sigma", "mu:sigma" = "mu:sigma", "sigma" = "sigma")
  if (is.null(which)) which <- names(p)
  
  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ## function for computing Hessian elements (expected or observed)
  hess <- if (expected) {
    function(par) switch(par,
      "mu"    = rep_len(-1/d$sigma^2, n),
      "sigma" = rep_len(-2 /d$sigma^2, n),
      rep.int(0, n))
  } else {
    function(par) switch(par,
      "mu"    = rep_len(-1/d$sigma^2, n),
      "sigma" = -3 * (x - d$mu)^2/(d$sigma^4) + 1/d$sigma^2,
      -2 * (x - d$mu)/d$sigma^3)
  }
  
  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    h <- setNames(hess(w), names(d))
  } else {
    h <- lapply(w, hess)
    h <- do.call("cbind", h)
    dimnames(h) <- list(names(d), w)
    if (!identical(w, which)) h <- h[, p[which], drop = FALSE]
    colnames(h) <- which
  }
  return(h)
}

## Poisson methods for score/hessian
score.Poisson <- function(d, x, which = "lambda", drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## only one parameter
  which <- match.arg(which, "lambda", several.ok = TRUE)

  ## compute score
  s <- x/d$lambda - 1
  if (!drop) s <- cbind("lambda" = s)
  return(s)
}

hessian.Poisson <- function(d, x, which = "lambda", drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## only one parameter
  which <- match.arg(which, "lambda", several.ok = TRUE)

  ## compute hessian
  h <- if (expected) rep_len(-1/d$lambda, n) else -x/d$lambda^2
  if (!drop) h <- cbind("lambda" = h)
  return(h)
}

## Bernoulli methods for score/hessian
score.Bernoulli <- function(d, x, which = "p", drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## only one parameter
  which <- match.arg(which, "p", several.ok = TRUE)

  ## compute score
  s <- (x - d$p)/(d$p * (1 - d$p))
  if (!drop) s <- cbind("p" = s)
  return(s)
}

hessian.Bernoulli <- function(d, x, which = "p", drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## only one parameter
  which <- match.arg(which, "p", several.ok = TRUE)

  ## compute hessian
  h <- if (expected) rep_len(-1/(d$p * (1 - d$p)), n) else -x/d$p^2 - (1 - x)/(1 - d$p)^2
  if (!drop) h <- cbind("p" = h)
  return(h)
}

## Binomial methods for score/hessian
score.Binomial <- function(d, x, which = "p", drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## only one parameter supported
  which <- match.arg(which, c("p", "size"), several.ok = TRUE)
  if (!identical(which, "p")) warning("only the scores with respect to 'p' are supported")

  ## compute score
  s <- (x - d$size * d$p)/(d$p * (1 - d$p))
  if (!drop) s <- cbind("p" = s)
  return(s)
}

hessian.Binomial <- function(d, x, which = "p", drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## only one parameter supported
  which <- match.arg(which, c("p", "size"), several.ok = TRUE)
  if (!identical(which, "p")) warning("only the scores with respect to 'p' are supported")

  ## compute hessian
  h <- if (expected) rep_len(-d$size/(d$p * (1 - d$p)), n) else -x/d$p^2 - (d$size - x)/(1 - d$p)^2
  if (!drop) h <- cbind("p" = h)
  return(h)
}


## for betareg
## BetaR methods for score/hessian

score.BetaR <- function(d, x, which = NULL, drop = TRUE, ...) {
  if (is.null(which)) which <- c("mu", "phi")
  s <- betareg:::sbetar(x, mu = d$mu, phi = d$phi, parameter = which, drop = FALSE)
  rownames(s) <- names(d)
  if (drop && NCOL(s) == 1L) s <- drop(s)
  return(s)
}

hessian.BetaR <- function(d, x, which = c("mu", "phi"), drop = TRUE, expected = TRUE, ...) {
  if (!identical(expected, TRUE)) stop("currently only the expected hessian is available")
  if (is.null(which)) which <- c("mu", "phi")
  h <- betareg:::hbetar(x, mu = d$mu, phi = d$phi, parameter = which, drop = FALSE) ## FIXME: support "phi:mu" in addition to "mu:phi"
  rownames(h) <- names(d)
  if (drop && NCOL(h) == 1L) h <- drop(h)
  return(-h) ## FIXME: erroneous sign in hbetar?
}


## for gamlss.dist
## GAMLSS methods for score/hessian

score.GAMLSS <- function(d, x, which = NULL, drop = TRUE, ...) {
  ## expand x if necessary
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)
  if (length(x) < n) x <- rep.int(x, n)

  ## get gamlss family object
  f <- gamlss.dist:::.gamlss_family(d)
  
  ## parameter names
  p <- names(f$parameters)
  if (is.null(which)) which <- p
  which <- match.arg(which, p, several.ok = TRUE)

  ## interface dl* functions from gamlss family
  dl <- c(mu = "dldm", sigma = "dldd", nu = "dldv", tau = "dldt")
  scr <- function(par) {
    dl_fun <- f[[dl[par]]]
    formals(dl_fun) <- c(formals(dl_fun), alist("..." = ))
    do.call(dl_fun, c(list(y = x), as.list(d), list(...)))
  }

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    s <- setNames(scr(which), names(d))
  } else {
    s <- lapply(which, scr)
    s <- do.call("cbind", s)
    dimnames(s) <- list(names(d), which)
  }
  return(s)
}

hessian.GAMLSS <- function(d, x, which = NULL, drop = TRUE, ...) {
  ## expand x if necessary
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)
  if (length(x) < n) x <- rep.int(x, n)

  ## get gamlss family object
  f <- gamlss.dist:::.gamlss_family(d)

  ## available and selected parameters/combinations and mappings for symmetries
  p <- names(f$parameters)
  pp <- outer(p, p, paste, sep = ":")
  diag(pp) <- p
  p <- setNames(
    c(diag(pp), pp[upper.tri(pp)], pp[upper.tri(pp)]),
    c(diag(pp), pp[upper.tri(pp)], pp[lower.tri(pp)])
  )[pp]
  if (is.null(which)) which <- names(p)
  
  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ## function for computing Hessian elements
  d2l <- c("mu" = "d2ldm2", "sigma" = "d2ldd2", "nu" = "d2ldv2", "tau" = "d2ldt2",
    "mu:sigma" = "d2ldmdd", "mu:nu" = "d2ldmdv", "mu:tau" = "d2ldmdt",
    "sigma:nu" = "d2ldddv", "sigma:tau" = "d2ldddt", "nu:tau" = "d2ldvdt")
  hess <- function(par) {
    d2l_fun <- f[[d2l[par]]]
    formals(d2l_fun) <- c(formals(d2l_fun), alist("..." = ))
    do.call(d2l_fun, c(list(y = x), as.list(d), list(...)))
  }
  
  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    h <- setNames(hess(w), names(d))
  } else {
    h <- lapply(w, hess)
    h <- do.call("cbind", h)
    dimnames(h) <- list(names(d), w)
    if (!identical(w, which)) h <- h[, p[which], drop = FALSE]
    colnames(h) <- which
  }
  return(h)
}
