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
distributions3_family <- function(distribution, links, score = TRUE, hessian = FALSE, update = FALSE, create_distribution = "structure",
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
    "log_likelihood"   = function(par, y, ...) sum(log_pdf(d3(par), y, elementwise = TRUE, ...)),
    "mu"       = function(par, ...) mean(d3(par), ...),
    "pdf"      = function(par, y, log = FALSE) pdf(d3(par), y, elementwise = TRUE, log = log),
    "cdf"      = function(par, y, ...) cdf(d3(par), y, elementwise = TRUE, ...),
    "random"   = function(par, n) random(d3(par), n),
    "quantile" = function(par, p) quantile(d3(par), p, elementwise = TRUE),
    "crps"     = function(par, y, ...) sum(scoringRules::crps(d3(par), y, elementwise = TRUE, ...)),
    "mean"     = function(par) mean(d3(par)),
    "variance" = function(par) variance(d3(par)),
    "skewness" = function(par) skewness(d3(par)),
    "kurtosis" = function(par) kurtosis(d3(par)),
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
  if (hessian && has_hessian) rval$hessian <- structure(lapply(nams, function(n) { ## FIXME: hessian(par, y, which, ...)
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
  if (update && has_score && has_hessian) rval$update <- function(par, y, eta, which) { ## FIXME: update(par, y, eta, which) -> list(eta, weights)
    lnk <- links[[which]]
    mu.eta <- lnk$mu.eta(eta)
    par <- d3(par)

    s <- score(par, y, which = which)
    score <- s * mu.eta
    score <- deriv_checks(score, is.weight = FALSE)

    h <- hessian(par, y, which = which)
    hessian <- h * mu.eta^2 + s * lnk$mu.eta2(eta)
    hessian <- deriv_checks(hessian, is.weight = TRUE)

    list(
      eta = eta + 1/hessian * score,
      weights = hessian
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
