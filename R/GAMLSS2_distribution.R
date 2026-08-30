## create GAMLSS2 (gamlss2 distribution family) distributions3 objects:
## a single class and corresponding methods
## encompassing all distributions using the workflow from
## the distributions3 package
GAMLSS2 <- function(family, parameters) {
  ## get family object
  f <- complete_family(family)
  ## set up distribution
  if(!is.list(parameters) && (length(f$names) > 1)) {
    if(is.null(dim(parameters)))
      parameters <- as.data.frame(t(parameters))
  }
  if(!inherits(parameters, "data.frame"))
    parameters <- as.data.frame(parameters)
  names(parameters) <- f$names
  class(parameters) <- c("GAMLSS2", "distribution")
  attr(parameters, "family") <- f
  return(parameters)
}

## S3 method for extracting fitted/predicted distributions3 objects
## associated methods are in gamlss.dist (as well as distributions3, topmodels, etc.)
prodist.gamlss2 <- function(object, ...) {
  ## extract fitted parameters
  d <- predict(object, ..., type = "parameter", drop = FALSE)

  ## set class
  if(is.null(object$family$create_distribution)) {
    class(d) <- c("GAMLSS2", "distribution")
  } else {
    d <- object$family$create_distribution(d)
  }

  ## include family information
  attr(d, "family") <- family(object)

  ## return distributions3 object
  return(d)
}

## S3 methods
format.GAMLSS2 <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
  class(x) <- c(paste("GAMLSS2", attr(x, "family")$family[1L]), "distribution")
  NextMethod()
}

print.GAMLSS2 <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
  class(x) <- c(paste("GAMLSS2", attr(x, "family")$family[1L]), "distribution")
  NextMethod()
}

mean.GAMLSS2 <- function(x, ...) {
  f <- complete_family(x)
  if(is.null(f$mean))
    stop(sprintf("the mean is not implemented for the %s family", attr(x, "family")[1L]))
  m <- f$mean(as.list(x))
  setNames(m, names(x))
}

variance.GAMLSS2 <- function(x, ...) {
  f <- complete_family(x)
  if(is.null(f$variance))
    stop(sprintf("the variance is not implemented for the %s family", attr(x, "family")[1L]))
  m <- f$variance(as.list(x))
  setNames(m, names(x))
}

skewness.GAMLSS2 <- function(x, ...) {
  f <- complete_family(x)
  if(is.null(f$skewness))
    stop(sprintf("the skewness is not implemented for the %s family", attr(x, "family")[1L]))
  m <- f$skewness(as.list(x))
  setNames(m, names(x))
}

kurtosis.GAMLSS2 <- function(x, ...) {
  f <- complete_family(x)
  if(is.null(f$kurtosis))
    stop(sprintf("the kurtosis is not implemented for the %s family", attr(x, "family")[1L]))
  m <- f$kurtosis(as.list(x))
  setNames(m, names(x))
}

pdf.GAMLSS2 <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  f <- complete_family(d)
  FUN <- function(at, d) { f$pdf(par = d, y = at, log = FALSE) }
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.GAMLSS2 <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  f <- complete_family(d)
  FUN <- function(at, d) { f$pdf(par = d, y = at, log = TRUE) }
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.GAMLSS2 <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  f <- complete_family(d)
  FUN <- function(at, d) { f$cdf(par = d, y = at) }
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.GAMLSS2 <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  f <- complete_family(x)
  FUN <- function(at, d) { f$quantile(d, at) }
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

random.GAMLSS2 <- function(x, n = 1L, drop = TRUE, ...) {
  f <- complete_family(x)
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) {
    ## apply_dpqr() handles replication for vectorized distributions and
    ## expects one draw per parameter row in that case. For a scalar
    ## distribution, however, it delegates all requested draws at once.
    f$random(d, if(length(d) == 1L) at else 1L)
  }
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

support.GAMLSS2 <- function(d, drop = TRUE, ...) {
  s <- quantile(d, probs = c(0, 1), elementwise = FALSE)
  distributions3::make_support(s[, 1L], s[, 2L], d, drop = drop)
}

is_discrete.GAMLSS2 <- function(d, ...) {
  f <- complete_family(d)
  if(is.null(f$type)) stop(sprintf("the type is not implemented for the %s family", attr(d, "family")[1L]))
  setNames(rep.int(tolower(f$type) == "discrete", length(d)), names(d))
}

is_continuous.GAMLSS2 <- function(d, ...) {
  f <- complete_family(d)
  if(is.null(f$type)) stop(sprintf("the type is not implemented for the %s family", attr(d, "family")[1L]))
  setNames(rep.int(tolower(f$type) == "continuous", length(d)), names(d))
}

