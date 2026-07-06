\name{deviance.gamlss2}
\alias{deviance.gamlss2}
\title{Deviance of GAMLSS Models}

\description{
  Extracts the deviance from fitted \code{gamlss2} objects. By default, the
  total deviance is returned. Optionally, observation-wise deviance contributions
  are returned.
}

\usage{
\method{deviance}{gamlss2}(object, \dots, newdata = NULL, sum = TRUE)
}

\arguments{
  \item{object}{An object of class \code{"gamlss2"}.}
  \item{\dots}{Optionally, further fitted \code{gamlss2} objects.}
  \item{newdata}{An optional data frame in which to evaluate the deviance. If
    omitted, the deviance is evaluated for the data used to fit the model.}
  \item{sum}{Should the observation-wise deviance contributions be summed? If
    \code{TRUE}, the default, the total deviance is returned. If \code{FALSE},
    the individual deviance contributions are returned.}
}

\details{
  The deviance of a fitted \code{gamlss2} model is defined as minus twice the
  fitted log-likelihood,
  \deqn{-2 \log L.}
  If \code{sum = TRUE}, one deviance value is returned for each supplied model.
  If more than one model is supplied, the deviances are ordered increasingly.

  If \code{sum = FALSE}, the function returns the observation-wise deviance
  contributions,
  \deqn{-2 \log f(y_i | \hat\theta_i),}
  where \eqn{f()} is the fitted probability density or probability mass function
  and \eqn{\hat\theta_i} denotes the fitted distribution parameters for
  observation \eqn{i}. These values sum to the total deviance and are useful, for
  example, in cross-validation.
}

\value{
  If \code{sum = TRUE}, a numeric vector with the total deviance of the supplied
  model or models. If \code{sum = FALSE}, a numeric vector of observation-wise
  deviance contributions for one model, or a matrix with one column per model if
  several models are supplied.
}

\seealso{
  \code{\link{gamlss2}}
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
## load the abdominal circumference data
data("abdom", package = "gamlss.data")

## specify a distributional model:
## location and scale vary with x,
## shape parameters remain constant
f <- y ~ s(x) | s(x) | 1 | 1

## estimate model
m1 <- gamlss2(f, data = abdom, family = BCT)

## total deviance
deviance(m1)

## observation-wise deviance contributions
di <- deviance(m1, sum = FALSE)

## the contributions sum to the total deviance
sum(di)
deviance(m1)

## compare two models
m2 <- gamlss2(y ~ s(x) | 1 | 1 | 1, data = abdom, family = BCT)
deviance(m1, m2)

## compare observation-wise deviance contributions
di <- deviance(m1, m2, sum = FALSE)
head(di)

## deviance on new data
nd <- abdom[1:20, ]
deviance(m1, newdata = nd)
deviance(m1, newdata = nd, sum = FALSE)
}

\keyword{models}
\keyword{regression}

