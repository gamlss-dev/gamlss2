\name{predict.gamlss2}

\alias{predict.gamlss2}

\title{Extracting Fitted or Predicted Parameters or Terms from gamlss2 Models}

\description{
Methods for \pkg{gamlss2} model objects that extract fitted values or predictions
for distribution parameters, additive predictors, individual terms, or the mean
of the response distribution.
}

\usage{
\method{predict}{gamlss2}(object, model = NULL, newdata = NULL,
  type = c("parameter", "link", "response", "terms"), terms = NULL,
  se.fit = FALSE, drop = TRUE, ...)
}
\arguments{
  \item{object}{A model object of class \code{\link{gamlss2}}.}
  \item{model}{Character. Which distribution parameter(s) should be predicted?
    This can be one or more of \code{"mu"}, \code{"sigma"}, etc. By default,
    predictions are returned for all model parts.}
  \item{newdata}{Data frame. Optionally, a new data frame in which to
    look for variables with which to predict. If omitted, the original
    observations are used.}
  \item{type}{Character. Which type of prediction should be computed?
    Use \code{"link"} for additive predictors on the link scale,
    \code{"parameter"} for distribution parameters on their natural scale,
    \code{"terms"} for individual additive terms, and \code{"response"} for
    the mean of the fitted response distribution.}
  \item{terms}{Character. Which of the terms in the additive predictor(s)
    should be included? By default all terms are included.}
  \item{se.fit}{Logical. Should standard errors for the predictions be included?
    Standard errors are computed by simulating from the approximate multivariate normal
    distribution of the maximum likelihood estimates. The number of simulations is controlled
    by the argument \code{R}, which defaults to \code{R = 200}, and can be passed via \code{...}.}
  \item{drop}{Logical. Should the predictions be simplified to a vector
    if possible (\code{TRUE}) or always returned as a data frame (\code{FALSE})?}
  \item{...}{Further control arguments. Supported aliases are \code{what} and
    \code{parameter} for \code{model}. Additional arguments are used for
    simulation-based summaries, including \code{FUN}, \code{R}, \code{seed},
    \code{burnin}, and \code{nogrep}.}
}

\details{
Prediction for a \code{\link{gamlss2}} object proceeds in stages. First, the
required model terms are evaluated for the original data or for \code{newdata}.
These terms are then added to form the full additive predictor for each
distribution parameter. Applying the inverse link function yields the predicted
distribution parameters. If requested, the mean of the corresponding response
distribution is computed in a final step.

If \code{se.fit = TRUE}, or if a summary function is supplied via \code{FUN},
predictions are based on simulation draws from the approximate multivariate
normal distribution of the fitted coefficients, unless the object already
contains stored draws in \code{samples}. For objects of class \code{"bamlss2"},
stored posterior draws are used directly.

For \code{type = "terms"}, the \code{terms} argument is matched against linear
and special-term labels. By default, partial matching is used; set
\code{nogrep = TRUE} in \code{...} to require exact matching.

See also \code{\link{prodist.gamlss2}} for setting up a full \pkg{distributions3} object
from which moments, probabilities, quantiles, or random numbers can be obtained.
}

\value{
If \code{drop = FALSE}, a data frame is returned. If \code{drop = TRUE}
 (the default), the result may be simplified to a numeric vector when possible.
 For \code{type = "terms"}, a data frame is returned with one column per selected
 term. For \code{se.fit = TRUE}, the returned object contains fitted values and
 corresponding simulation-based standard errors.
}

\seealso{
\code{\link[stats]{predict}}, \code{\link{prodist.gamlss2}}, \code{\link{quantile.gamlss2}}
}

\examples{
## fit heteroscedastic normal GAMLSS model
## stopping distance (ft) explained by speed (mph)
data("cars", package = "datasets")
m <- gamlss2(dist ~ s(speed) | s(speed), data = cars, family = NO)

## new data for predictions
nd <- data.frame(speed = c(10, 20, 30))

## default: model parameter(s) for all model parts
predict(m, newdata = nd)

## additive predictors on the link scale
predict(m, newdata = nd, type = "link")

## mean of the response distribution
predict(m, newdata = nd, type = "response")

## model parameter(s)
predict(m, newdata = nd, model = "sigma")
predict(m, newdata = nd, model = "sigma", drop = FALSE)

## individual terms in additive predictor(s)
predict(m, newdata = nd, type = "terms", model = "sigma")
predict(m, newdata = nd, type = "terms", model = "sigma", terms = "s(speed)")

## predict quantiles
quantile(m, newdata = nd, probs = c(0.1, 0.5, 0.9))

## standard errors
predict(m, newdata = nd, se.fit = TRUE, R = 200)
}

\keyword{regression}
