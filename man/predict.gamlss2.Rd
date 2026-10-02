\name{predict.gamlss2}

\alias{predict.gamlss2}

\title{Predictions from gamlss2 Models}

\description{
Predict response distributions, parameters, additive predictors, terms, or
response means from a \pkg{gamlss2} model.
}

\usage{
\method{predict}{gamlss2}(object, parameter = NULL, newdata = NULL,
  type = c("distribution", "parameter", "link", "response", "terms"),
  terms = NULL, se.fit = FALSE, drop = TRUE, ...,
  level = NULL, interval = c("none", "wald"), unconditional = FALSE)
}

\arguments{
  \item{object}{A \code{\link{gamlss2}} model.}
  \item{parameter}{Distribution parameter name(s) or index(es), such as
    \code{"mu"} or \code{"sigma"}; defaults to all parameters.}
  \item{newdata}{Data frame for predictions; defaults to the fitted data.}
  \item{type}{\code{"distribution"} (default) returns \pkg{distributions3}
    objects; \code{"parameter"} returns parameters on their natural scale;
    \code{"link"} returns additive predictors; \code{"response"} returns
    response means; \code{"terms"} returns individual additive terms.
    Specifying \code{parameter} without \code{type} selects
    \code{"parameter"}.}
  \item{terms}{Additive terms to include; defaults to all. Names are partially
    matched unless \code{nogrep = TRUE} is passed in \code{...}.}
  \item{se.fit}{For ordinary predictions, \code{TRUE} computes standard errors
    from coefficient draws; \code{R} in \code{...} sets the number of draws
    (default 200).}
  \item{drop}{Simplify the result when possible.}
  \item{level}{Numeric confidence level(s) between zero and one; supplying
    levels selects pointwise Wald intervals.}
  \item{interval}{\code{"none"} (default) returns ordinary predictions;
    \code{"wald"} returns analytic pointwise confidence intervals, at the 95\%
    confidence level by default.}
  \item{unconditional}{If \code{TRUE}, include approximate smoothing-parameter
    uncertainty in Wald intervals and coefficient draws.}
  \item{...}{Further controls, including \code{FUN}, \code{R}, \code{seed},
    and \code{burnin} for simulation.}
}

\details{
Coefficient draws stored in the fit are used when available; otherwise, draws
use the full joint coefficient covariance. They approximate the coefficient
sampling distribution and may give unreliable tails after nonlinear
transformations. Simulation stops if a marginal 95\% Gaussian interval for a
log smoothing parameter spans a factor of 100 on the smoothing-parameter scale.

Wald intervals are approximate pointwise \emph{confidence} intervals for fitted
values, not \emph{prediction} intervals for future observations. They are
computed without simulation and may be asymmetric on transformed scales.
They condition on fitted smoothing parameters by default and exclude
model-selection uncertainty; see \code{\link{vcov.gamlss2}}.

Wald intervals support maximum-likelihood fits with linear terms and
coefficient-linear \pkg{mgcv} smooths, but not nonlinear custom special terms.
They cannot be combined with \code{FUN}, \code{R}, \code{seed}, or \code{burnin}.
An indefinite joint information matrix or non-estimable coefficient causes an
error rather than a working-information approximation.
}

\value{
Ordinary predictions return a distribution object or values in the selected
form. Term predictions are data frames or, for multiple parameters, lists of
matrices. With \code{se.fit = TRUE}, fitted values and standard errors are
returned. With \code{FUN = identity}, simulated distribution predictions form
a list of distribution objects.

Wald predictions return \code{fit}, \code{se.fit}, \code{lower}, and \code{upper}.
For multiple levels, \code{lower} and \code{upper} are lists named by confidence
percentage (for example, \code{"80\%"} and \code{"95\%"}). \code{fit} and
each bound retain the structure selected by \code{parameter}, \code{terms}, and
\code{drop}; \code{se.fit} is shared across levels.
}

\seealso{
\code{\link[stats]{predict}}, \code{\link{marginal_predict}},
\code{\link{prodist.gamlss2}}, \code{\link{quantile.gamlss2}}
}

\examples{
## fit heteroscedastic normal GAMLSS model
## stopping distance (ft) explained by speed (mph)
data("cars", package = "datasets")
m <- gamlss2(dist ~ s(speed) | s(speed), data = cars, family = NO)

## new data for predictions
nd <- data.frame(speed = c(10, 20, 30))

## default: a "distribution" object
d <- predict(m, newdata = nd)
print(d)

## can be used with distributions3 methods
pdf(d, 20)
cdf(d, 20)
quantile(d, 0.8)
mean(d)
variance(d)
random(d, 3)

## additive predictors on the link scale
predict(m, newdata = nd, type = "link")

## mean of the response distribution
predict(m, newdata = nd, type = "response")

## same as
mean(d)

## all model parameters
predict(m, newdata = nd, type = "parameter")

## parameter selects parameter predictions if type is omitted
predict(m, newdata = nd, parameter = "sigma")
predict(m, newdata = nd, parameter = "sigma", drop = FALSE)

## individual terms in additive predictor(s)
predict(m, newdata = nd, type = "terms", parameter = "sigma")
predict(m, newdata = nd, type = "terms", parameter = "sigma", terms = "s(speed)")

## predict quantiles
quantile(m, newdata = nd, probs = c(0.1, 0.5, 0.9))

## same as
quantile(d, probs = c(0.1, 0.5, 0.9), elementwise = FALSE)

## standard errors
predict(m, type = "parameter", newdata = nd, se.fit = TRUE, R = 200)

## fast pointwise confidence bands, without sampling
bands <- predict(m, type = "parameter", newdata = nd, level = c(0.8, 0.95))
bands$lower[["80\%"]]
bands$upper[["95\%"]]

## prediction uncertainty using multiple draws
d <- predict(m, newdata = cars, R = 1000, FUN = identity)

## predict 1000 quantiles
q <- sapply(d, quantile, probs = 0.95)

i95 <- function(x) {
  setNames(quantile(x, c(0.025, 0.5, 0.975), names = FALSE),
    c("lower", "fit", "upper"))
}

## quantiles of quantile samples
q <- t(apply(q, 1, i95))

## plot quantile with uncertainty bands
plot(dist ~ speed, data = cars,
  ylim = range(cars$dist, q))
matlines(cars$speed, q,
  lty = c(2, 1, 2), col = 4, lwd = 2)
}

\keyword{regression}
