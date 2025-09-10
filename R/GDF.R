## create GDF (gamlss2 distribution family) distributions3 objects:
## a single class and corresponding methods
## encompassing all distributions using the workflow from
## the distributions3 package
GDF <- function(family, parameters) {
  stopifnot(requireNamespace("distributions3"))
  ## get family object
  f <- .GDF_family(family)
  ## set up distribution
  if(!is.list(parameters) && (length(f$names) > 1)) {
    if(is.null(dim(parameters)))
      parameters <- as.data.frame(t(parameters))
  }
  if(!inherits(parameters, "data.frame"))
    parameters <- as.data.frame(parameters)
  names(parameters) <- f$names
  class(parameters) <- c("GDF", "distribution")
  attr(parameters, "family") <- f
  return(parameters)
}

## S3 method for extracting fitted/predicted distributions3 objects
## associated methods are in gamlss.dist (as well as distributions3, topmodels, etc.)
prodist.gamlss2 <- function(object, ...) {
  stopifnot(requireNamespace("distributions3"))

  ## extract fitted parameters
  d <- predict(object, ..., type = "parameter", drop = FALSE)

  ## set class to general GAMLSS or GDF distribution (distributions3 object)
  if(family(object)$family[1L] %in% getNamespaceExports("gamlss.dist") && FALSE) {
    class(d) <- c("GAMLSS", "distribution")
  } else {
    class(d) <- c("GDF", "distribution")
  }

  ## include family information
  attr(d, "family") <- family(object)

  ## return distributions3 object
  return(d)
}

## auxiliary functions for getting the family.
.GDF_family <- function(family) {
  if(inherits(family, "GDF")) {
    family <- attr(family, "family")
  }
  if(is.character(family)) {
    family <- get(family)
  }
  if(is.function(family)) {
    family <- family()
  }
  if(inherits(family, "gamlss.family")) {
    family <- tF(family)
  }
  if(is.null(family$type)) {
    warning(paste0('no type specified in family ', family$family[1L],
      ', setting type = "Continuous"!'))
    family$type <- "Continuous"
  }
  cf <- class(family)
  elmts <- c(
    "family", "names", "pdf", "cdf", "quantile", "random",
    "mean", "mode", "variance", "skewness", "kurtosis",
    "type"
  )
  family <- family[elmts]
  family <- family[!unlist(lapply(family, is.null))]
  class(family) <- cf
  return(family)
}

## S3 methods
format.GDF <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
  stopifnot(requireNamespace("distributions3"))
  class(x) <- c(paste("GDF", attr(x, "family")$family[1L]), "distribution")
  NextMethod()
}

print.GDF <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
  stopifnot(requireNamespace("distributions3"))
  class(x) <- c(paste("GDF", attr(x, "family")$family[1L]), "distribution")
  NextMethod()
}

mean.GDF <- function(x, ...) {
  f <- .GDF_family(x)
  if(is.null(f$mean))
    stop(sprintf("the mean is not implemented for the %s family", attr(x, "family")[1L]))
  m <- f$mean(as.list(x))
  setNames(m, names(x))
}

variance.GDF <- function(x, ...) {
  f <- .GDF_family(x)
  if(is.null(f$variance))
    stop(sprintf("the variance is not implemented for the %s family", attr(x, "family")[1L]))
  m <- f$variance(as.list(x))
  setNames(m, names(x))
}

skewness.GDF <- function(x, ...) {
  f <- .GDF_family(x)
  if(is.null(f$skewness))
    stop(sprintf("the skewness is not implemented for the %s family", attr(x, "family")[1L]))
  m <- f$skewness(as.list(x))
  setNames(m, names(x))
}

kurtosis.GDF <- function(x, ...) {
  f <- .GDF_family(x)
  if(is.null(f$kurtosis))
    stop(sprintf("the kurtosis is not implemented for the %s family", attr(x, "family")[1L]))
  m <- f$kurtosis(as.list(x))
  setNames(m, names(x))
}

pdf.GDF <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  f <- .GDF_family(d)
  FUN <- function(at, d) { f$pdf(at, d, log = FALSE) }
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.GDF <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  f <- .GDF_family(d)
  FUN <- function(at, d) { f$pdf(at, d, log = TRUE) }
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.GDF <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  f <- .GDF_family(d)
  FUN <- function(at, d) { f$cdf(at, d) }
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.GDF <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  f <- .GDF_family(x)
  FUN <- function(at, d) { f$quantile(at, d) }
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

random.GDF <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  f <- .GDF_family(x)
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) { f$random(at, d) }
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

support.GDF <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  s <- quantile(d, probs = c(0, 1), elementwise = FALSE)
  distributions3::make_support(s[, 1L], s[, 2L], d, drop = drop)
}

is_discrete.GDF <- function(d, ...) {
  f <- .GDF_family(d)
  if(is.null(f$type)) stop(sprintf("the type is not implemented for the %s family", attr(d, "family")[1L]))
  setNames(rep.int(tolower(f$type) == "discrete", length(d)), names(d))
}

is_continuous.GDF <- function(d, ...) {
  f <- .GDF_family(d)
  if(is.null(f$type)) stop(sprintf("the type is not implemented for the %s family", attr(d, "family")[1L]))
  setNames(rep.int(tolower(f$type) == "continuous", length(d)), names(d))
}

