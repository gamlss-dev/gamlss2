\name{logLik.gamlss2}
\alias{logLik.gamlss2}
\title{Log-Likelihood of GAMLSS Models}

\description{
  Extracts the fitted log-likelihood from \code{gamlss2} objects. Optionally,
  the log-likelihood can be evaluated on new data. Several fitted models can be
  supplied for comparison.
}

\usage{
\method{logLik}{gamlss2}(object, \dots, newdata = NULL)
}

\arguments{
  \item{object}{An object of class \code{"gamlss2"}.}
  \item{\dots}{Optionally, further fitted \code{gamlss2} objects.}
  \item{newdata}{An optional data frame in which to evaluate the
    log-likelihood. If omitted, the log-likelihood stored in the fitted model
    object is returned. If supplied, \code{newdata} must contain the response
    variable.}
}

\details{
  For the original model data, \code{logLik()} returns the log-likelihood value
  stored in the fitted \code{gamlss2} object. If \code{newdata} is supplied, the
  distribution parameters are first predicted for the observations in
  \code{newdata}; the response is then extracted from \code{newdata} and the
  log-likelihood is recomputed using the model family.

  The returned object follows the usual \code{\link[stats]{logLik}} convention:
  it has class \code{"logLik"} and carries attributes \code{"df"} and
  \code{"nobs"}. The degrees of freedom are taken from the fitted model object,
  while the number of observations is either the number of observations used for
  fitting or the number of rows evaluated in \code{newdata}.

  If more than one fitted model is supplied, a named list of \code{"logLik"}
  objects is returned. The list is ordered by decreasing log-likelihood.
}

\value{
  For a single model, an object of class \code{"logLik"}. For several supplied
  models, a named list of \code{"logLik"} objects, ordered by decreasing
  log-likelihood.
}

\seealso{
  \code{\link[stats]{logLik}}, \code{\link{deviance.gamlss2}},
  \code{\link{gamlss2}}
}

\examples{
## fit two normal GAMLSS models
## stopping distance (ft) explained by speed (mph)
data("cars", package = "datasets")

m1 <- gamlss2(dist ~ speed | 1, data = cars, family = NO)
m2 <- gamlss2(dist ~ speed | speed, data = cars, family = NO)

## fitted log-likelihood
logLik(m1)

## attributes follow the usual logLik convention
ll <- logLik(m1)
attr(ll, "df")
attr(ll, "nobs")

## compare several fitted models
logLik(m1, m2)

## log-likelihood evaluated on new data
nd <- cars[1:10, ]
logLik(m2, newdata = nd)

## relation to the deviance
-2 * as.numeric(logLik(m2))
deviance(m2)
}

\keyword{models}
\keyword{regression}

