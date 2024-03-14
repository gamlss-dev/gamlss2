\name{predict.gamlss2}

\alias{predict.gamlss2}

\title{Extracting Fitted or Predicted Parameters or Terms from gamlss2 Models}

\description{
Methods for \pkg{gamlss2} model objects for extracting fitted (in-sample) or
predicted (out-of-sample) parameters, terms, etc.
}

\usage{
\method{predict}{gamlss2}(object, model = NULL, newdata = NULL,
  type = c("parameter", "link", "response", "terms"), terms = NULL,
  se.fit = FALSE, drop = TRUE, ...)
}
\arguments{
  \item{object}{model object of class \code{\link{gamlss2}}.}
  \item{model}{character. Which model part(s) should be predicted?
    Can be one or more of \code{"mu"}, \code{"sigma"}, etc. By default
    all model parts are included.}
  \item{newdata}{data.frame. Optionally, a new data frame in which to
    look for variables with which to predict. If omitted, the original
    observations are used.}
  \item{type}{character. Which type of prediction should be computed?
    Can be the full additive predictor(s) (\code{"link"}, before applying the link function(s)),
    the corresponding parameter (\code{"parameter"}, after applying the link function(s)),
    the individual terms of the additive predictor(s) (\code{"terms"}),
    or the corresponding mean of the response distribution (\code{"response"}).}
  \item{terms}{character. Which of the terms in the additive predictor(s)
    should be included? By default all terms are included.}
  \item{se.fit}{logical. Should standard errors for the predictions be included?
    \emph{(not implemented yet)}.}
  \item{drop}{logical. Should the predictions be simplified to a vector
    if possible (\code{TRUE}) or always returned as a data.frame (\code{FALSE})?}
  \item{...}{currently only used for catching \code{what} as an alias for \code{model}.}
}

\details{
Predictions for \code{\link{gamlss2}} model objects are obtained in the following steps:
First, the original data is extracted or some \code{newdata} is set up.
Second, all of the terms in the additive predictors of all model parameters
(\code{"mu"}, \code{"sigma"}, ...) are computed.
Third, the full additive predictor(s) are obtained by adding up all individual terms.
Fourth, the parameter(s) are obtained from the additive predictor(s) by applying the inverse
link function(s).
In a final step, the mean of the associated probability distribution can be computed.

See also \code{\link{prodist.gamlss2}} for setting up a full \pkg{distributions3} object
from which moments, probabilities, quantiles, or random numbers can be obtained.
}

\value{
If \code{drop = FALSE} a data.frame. If \code{drop = TRUE} (the default),
the data.frame might be simplified to a numeric vector, if possible.
}

\seealso{
\code{\link[stats]{predict}}, \code{\link{prodist.gamlss2}}
}

\examples{
## fit heteroscedastic normal GAMLSS model
## stopping distance (ft) explained by speed (mph)
data("cars", package = "datasets")
m <- gamlss2(dist ~ s(speed) | s(speed), data = cars, family = NO)

## new data for predictions
nd <- data.frame(speed = c(10, 20, 30))

## default: additive predictors (on link scale) for all model parameters
predict(m, newdata = nd)

## mean of the response distribution
predict(m, newdata = nd, type = "response")

## model parameter(s)
predict(m, newdata = nd)
predict(m, newdata = nd, model = "sigma")
predict(m, newdata = nd, model = "sigma", drop = FALSE)

## individual terms in additive predictor(s)
predict(m, newdata = nd, type = "terms", model = "sigma")
predict(m, newdata = nd, type = "terms", model = "sigma", terms = "s(speed)")
}

\keyword{regression}
