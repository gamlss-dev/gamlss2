## auxiliary function for checking whether dedicated S3 methods are available
hasS3method <- function (method, classes) {
  any(sapply(classes, function(cls) {
    tryCatch(is.function(getS3method(method, class = cls)), error = function(e) FALSE)
  }))
}

## S3 method for generating a gamlss2.family from a distribution object (from distributions3)
family.distribution <- function(object, links, score = TRUE, hessian = FALSE, update = FALSE, expected = FALSE, ...) {
  ## distributions3 class
  d <- setdiff(class(object), "distribution")
  create_distribution <- function(par) structure(par, class = class(object))

  ## score and hessian methods available?
  has_score <- hasS3method("score", d)
  has_hessian <- hasS3method("hessian", d)

  ## process links
  if(!is.list(links))
    links <- as.list(links)
  if(is.null(names(links))) {
    names(links) <- names(as.list(object))
  }
  
  distributions3_family(d[1L],
    links = links,
    score = score && has_score,
    hessian = hessian && has_hessian,
    update = update && has_score && has_hessian,
    expected = expected,
    create_distribution = create_distribution,
    ...)
}

## S3 methods for generating a gamlss2.family from gamlss.dist objects
## (either a gamlss.family or a GAMLSS distribution)
family.GAMLSS <- function(object, ...) {
  .gamlss_family <- function (x) get(attr(object, "family")[1L], asNamespace("gamlss.dist"))()
  object <- .gamlss_family(object)
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
distributions3_family <- function(distribution, links, score = TRUE, hessian = FALSE, update = FALSE, expected = FALSE, create_distribution = "structure",
  type = NULL, valid.response = NULL, initialize = NULL, ...) {
  ## distribution must be a character,
  ## corresponding to a distribution generator
  stopifnot(is.character(distribution))
  args <- argsAnywhere(distribution)

  ## If args is a list (found the same function in multiple packages) but
  ## one is a function from the distributions3 package, pick the distributions3
  ## version.
  if (is.list(args) && length(grep("^package:distributions3$", names(args), value = TRUE)) == 1L) {
      args <- args[[grep("^package:distributions3$", names(args), value = TRUE)[[1]]]]
  }
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
    "log_likelihood"   = function(par, y, ...) sum(log_pdf(d3(par), y, ...)),
    "mu"       = function(par, ...) mean(d3(par), ...),
    "pdf"      = function(par, y, log = FALSE) {
      if(log) log_pdf(d3(par), y) else pdf(d3(par), y)
    },
    "cdf"      = function(par, y, ...) cdf(d3(par), y, ...),
    "random"   = function(par, n) random(d3(par), n),
    "quantile" = function(par, p) quantile(d3(par), p),
    "crps"     = function(par, y, ...) sum(scoringRules::crps(d3(par), y, ...)),
    "mean"     = function(par) mean(d3(par)),
    "variance" = function(par) variance(d3(par)),
    "skewness" = function(par) skewness(d3(par)),
    "kurtosis" = function(par) kurtosis(d3(par)),
    "create_distribution" = d3
  )

  ## Prefer the native distribution support when its class provides one.
  if(hasS3method("support", distribution))
    rval$support <- function(par, ...) support(d3(par), drop = FALSE, ...)

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
  if (hessian && has_hessian) rval$hessian <- structure(lapply(nams, function(n) {
    if (expected) function(par, y, ...) {
      lnk <- links[[n]]
      eta <- lnk$linkfun(par[[n]])
      h <- hessian(d3(par), y, which = n, expected = TRUE, ...)
      h * lnk$mu.eta(eta)^2
    } else function(par, y, ...) {
      lnk <- links[[n]]
      eta <- lnk$linkfun(par[[n]])
      par <- d3(par)
      s <- score(par, y, which = n, ...)
      h <- hessian(par, y, which = n, ...)
      h * lnk$mu.eta(eta)^2 + s * lnk$mu.eta2(eta)
    }
  }), names = nams)

  ## add Newton/Fisher update function if desired and available
  if (update && has_score && has_hessian && expected) rval$update <- function(par, y, eta, which) { ## FIXME: update(par, y, eta, which) -> list(eta, weights)
    lnk <- links[[which]]
    mu.eta <- lnk$mu.eta(eta)
    par <- d3(par)

    s <- score(par, y, which = which)
    score <- s * mu.eta
    score <- deriv_checks(score, is.weight = FALSE)

    h <- hessian(par, y, which = which, expected = TRUE)
    hessian <- h * mu.eta^2
    hessian <- deriv_checks(hessian, is.weight = TRUE)

    list(
      eta = eta + 1/hessian * score,
      weights = hessian
    )
  }

  if (update && has_score && has_hessian && !expected) rval$update <- function(par, y, eta, which) {
    lnk <- links[[which]]
    mu.eta <- lnk$mu.eta(eta)
    par <- d3(par)
    s <- score(par, y, which = which)
    score <- deriv_checks(s * mu.eta, is.weight = FALSE)
    h <- hessian(par, y, which = which)
    hessian <- deriv_checks(h * mu.eta^2 + s * lnk$mu.eta2(eta), is.weight = TRUE)
    list(eta = eta + 1 / hessian * score, weights = hessian)
  }
  ## family elements that are difficult to infer automatically
  rval$valid.response <- if (is.null(valid.response)) function(x) TRUE else valid.response
  if (!is.null(type)) rval$type <- match.arg(type, c("continuous", "discrete")) ## FIXME: further types allowed?
  if (!is.null(initialize) && is.list(initialize) && all(names(initialize) %in% nams)) rval$initialize <- initialize

  ## return classed list
  class(rval) <- "gamlss2.family"
  return(rval)
}

##-----------------------------------------------------------------------------

family.Normal <- function(object, links = NULL, initialize = NULL, ...) {
  if (is.null(links)) links <- c(mu = "identity", sigma = "log")
  if (is.null(initialize)) initialize <- list(
    mu = function(y, ...) (y + mean(y))/2,
    sigma = function(y, ...) 0 * y + sd(y))
  if (isFALSE(initialize)) initialize <- NULL
  NextMethod(object, links = links, initialize = initialize, ...)
}

family.Poisson <- function(object, links = NULL, initialize = NULL, ...) {
  if (is.null(links)) links <- c(lambda = "log")
  if (is.null(initialize)) initialize <- list(lambda = function(y, ...) (y + mean(y))/2)
  if (isFALSE(initialize)) initialize <- NULL
  NextMethod(object, links = links, initialize = initialize, ...)
}

family.Bernoulli <- function(object, links = NULL, initialize = NULL, ...) {
  if (is.null(links)) links <- c(p = "logit")
  if (is.null(initialize)) initialize <- list(p = function(y, ...) (y + 0.5)/2)
  if (isFALSE(initialize)) initialize <- NULL
  NextMethod(object, links = links, initialize = initialize, type = "discrete",
    valid.response = function(x) all(is.finite(x) & x %in% c(0, 1)), ...)
}

family.NegativeBinomial <- function(object, links = NULL, initialize = NULL, ...) {
  nams <- names(as.list(object))
  if (!is.null(names(links))) nams <- names(links)
  nams <- if ("mu" %in% nams) c("mu", "size") else c("p", "size")
  if (is.null(links)) links <- if (nams[1L] == "mu") c(mu = "log", size = "log") else c(p = "logit", size = "log")
  if (is.null(initialize)) {
    ## Estimate the marginal size with a finite near-Poisson bound.
    initial <- function(y) {
      if (!any(y > 0)) return(100)
      mu <- max(mean(y), 1e-8)
      exp(optimize(function(eta) {
        -sum(stats::dnbinom(y, size = exp(eta), mu = mu, log = TRUE))
      }, interval = log(c(0.01, 100)), tol = 1e-6)$minimum)
    }
    initialize <- list(
      mu = function(y, ...) 0 * y + max(mean(y), 1e-8),
      p = function(y, ...) {
        size <- initial(y)
        0 * y + size/(size + max(mean(y), 1e-8))
      },
      size = function(y, ...) 0 * y + initial(y))
    initialize <- initialize[nams]
  }
  if (isFALSE(initialize)) initialize <- NULL
  rval <- NextMethod(object, links = links, initialize = initialize, type = "discrete",
    valid.response = function(x) all(is.finite(x) & x >= 0 & x == floor(x)), ...)
  ## Model the mean or probability before the size parameter.
  rval$names <- intersect(nams, rval$names)
  rval$links <- rval$links[rval$names]
  return(rval)
}

family.Gamma <- function(object, links = NULL, initialize = NULL, ...) {
  if (is.null(links)) links <- c(shape = "log", rate = "log")
  if (is.null(initialize)) {
    ## Log-moment shape estimate followed by two Newton corrections.
    initial <- function(y) {
      if (all(y == y[1L])) return(1)
      z <- log(y)
      d <- max(log(mean(exp(z - max(z)))) + max(z) - mean(z), 1e-8)
      shape <- (3 - d + sqrt((d - 3)^2 + 24 * d))/(12 * d)
      for (i in 1:2) {
        candidate <- shape - (log(shape) - digamma(shape) - d)/(1/shape - trigamma(shape))
        if (is.finite(candidate) && candidate > 0) shape <- candidate
      }
      return(shape)
    }
    initialize <- list(
      shape = function(y, ...) 0 * y + initial(y),
      rate = function(y, ...) 0 * y + initial(y)/mean(y))
  }
  if (isFALSE(initialize)) initialize <- NULL
  NextMethod(object, links = links, initialize = initialize, type = "continuous",
    valid.response = function(x) all(is.finite(x) & x > 0), ...)
}

family.LogNormal <- function(object, links = NULL, initialize = NULL, ...) {
  if (is.null(links)) links <- c(log_mu = "identity", log_sigma = "log")
  if (is.null(initialize)) initialize <- list(
    log_mu = function(y, ...) 0 * y + mean(log(y)),
    log_sigma = function(y, ...) 0 * y + max(sqrt(mean((log(y) - mean(log(y)))^2)), .Machine$double.eps^0.5))
  if (isFALSE(initialize)) initialize <- NULL
  NextMethod(object, links = links, initialize = initialize, type = "continuous",
    valid.response = function(x) all(is.finite(x) & x > 0), ...)
}

family.Beta <- function(object, links = NULL, initialize = NULL, ...) {
  if (is.null(links)) links <- c(alpha = "log", beta = "log")
  if (is.null(initialize)) {
    ## Moment estimates followed by five joint, damped Newton updates.
    initial <- function(y) {
      mu <- mean(y)
      v <- var(y)
      if (!is.finite(v) || v <= 0) v <- mu * (1 - mu)/3
      phi <- min(1e6, max(0.01, mu * (1 - mu)/v - 1))
      par <- c(mu, 1 - mu) * phi
      logs <- c(mean(log(y)), mean(log1p(-y)))
      ll <- function(par) lgamma(sum(par)) - sum(lgamma(par)) + sum((par - 1) * logs)
      value <- ll(par)
      for (i in 1:5) {
        s <- digamma(sum(par)) - digamma(par) + logs
        h <- diag(trigamma(par)) - trigamma(sum(par))
        delta <- tryCatch(drop(solve(h, s)), error = function(e) c(NA_real_, NA_real_))
        if (any(!is.finite(delta))) break
        accepted <- FALSE
        for (j in 0:19) {
          candidate <- par + delta/2^j
          if (all(is.finite(candidate)) && all(candidate > 0)) {
            next.value <- ll(candidate)
            if (is.finite(next.value) && next.value >= value - 1e-12) {
              par <- candidate
              value <- next.value
              accepted <- TRUE
              break
            }
          }
        }
        if (!accepted) break
      }
      return(par)
    }
    initialize <- list(
      alpha = function(y, ...) 0 * y + initial(y)[[1L]],
      beta = function(y, ...) 0 * y + initial(y)[[2L]])
  }
  if (isFALSE(initialize)) initialize <- NULL
  NextMethod(object, links = links, initialize = initialize, type = "continuous",
    valid.response = function(x) all(is.finite(x) & x > 0 & x < 1), ...)
}

family.Weibull <- function(object, links = NULL, initialize = NULL, ...) {
  if (is.null(links)) links <- c(shape = "log", scale = "log")
  if (is.null(initialize)) {
    ## Refine the log-moment shape and obtain its profile-likelihood scale.
    initial <- function(y) {
      z <- log(y)
      mu <- mean(z)
      v <- var(z)
      if (!is.finite(v) || v <= 0) v <- pi^2/6
      shape <- min(1e4, max(1e-4, pi/sqrt(6 * v)))
      z <- z - mu
      for (i in 1:3) {
        w <- shape * z
        w <- exp(w - max(w))
        w <- w/sum(w)
        m <- sum(w * z)
        s <- 1/shape - m
        h <- -1/shape^2 - sum(w * (z - m)^2)
        candidate <- shape - s/h
        if (is.finite(candidate) && candidate > 0) shape <- candidate
      }
      w <- shape * z
      scale <- exp(mu + (max(w) + log(mean(exp(w - max(w)))))/shape)
      return(c(shape, scale))
    }
    initialize <- list(
      shape = function(y, ...) 0 * y + initial(y)[[1L]],
      scale = function(y, ...) 0 * y + initial(y)[[2L]])
  }
  if (isFALSE(initialize)) initialize <- NULL
  NextMethod(object, links = links, initialize = initialize, type = "continuous",
    valid.response = function(x) all(is.finite(x) & x > 0), ...)
}

family.Binomial <- function(object, links = NULL, initialize = NULL, ...) {
  size <- object$size
  stopifnot(length(size) > 0L, all(is.finite(size) & size >= 0 & size == floor(size)))
  if (all(size == size[1L])) size <- size[1L]
  if (is.null(links)) links <- c(p = "logit")
  stopifnot(length(links) == 1L, is.null(names(links)) || identical(names(links), "p"))
  names(links) <- "p"
  if (is.null(initialize)) initialize <- list(p = function(y, ...) (y + 0.5)/(size + 1))
  if (isFALSE(initialize)) initialize <- NULL
  ## The number of trials is known and is not modeled.
  create_distribution <- function(par) distributions3::Binomial(size = size, p = par$p)
  rval <- distributions3_family("Binomial", links = links, initialize = initialize, type = "discrete",
    valid.response = function(x) {
      (length(size) == 1L || length(size) == length(x)) &&
        all(is.finite(x) & x >= 0 & x <= size & x == floor(x))
    }, create_distribution = create_distribution, ...)
  ## Avoid the legacy binomial response conversion in gamlss2().
  rval$family <- "distributions3::Binomial"
  return(rval)
}

family.Exponential <- function(object, links = NULL, initialize = NULL, ...) {
  if (is.null(links)) links <- c(rate = "log")
  if (is.null(initialize)) initialize <- list(
    rate = function(y, ...) {
      mu <- mean(y)
      if (mu == 0) mu <- 1e-8
      0 * y + 1/max(mu, .Machine$double.xmin)
    })
  if (isFALSE(initialize)) initialize <- NULL
  NextMethod(object, links = links, initialize = initialize, type = "continuous",
    valid.response = function(x) all(is.finite(x) & x >= 0), ...)
}

family.Logistic <- function(object, links = NULL, initialize = NULL, ...) {
  if (is.null(links)) links <- c(location = "identity", scale = "log")
  if (is.null(initialize)) initialize <- list(
    location = function(y, ...) (y + mean(y))/2,
    scale = function(y, ...) 0 * y + max(sd(y) * sqrt(3)/pi, .Machine$double.eps^0.5, na.rm = TRUE))
  if (isFALSE(initialize)) initialize <- NULL
  NextMethod(object, links = links, initialize = initialize, type = "continuous",
    valid.response = function(x) all(is.finite(x)), ...)
}

family.Cauchy <- function(object, links = NULL, initialize = NULL, ...) {
  if (is.null(links)) links <- c(location = "identity", scale = "log")
  if (is.null(initialize)) initialize <- list(
    location = function(y, ...) 0 * y + median(y),
    scale = function(y, ...) 0 * y + max(diff(quantile(y, c(0.25, 0.75)))/2, .Machine$double.eps^0.5))
  if (isFALSE(initialize)) initialize <- NULL
  NextMethod(object, links = links, initialize = initialize, type = "continuous",
    valid.response = function(x) all(is.finite(x)), ...)
}

