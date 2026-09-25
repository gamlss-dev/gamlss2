\name{predict.gamlss2}

\alias{predict.gamlss2}

\title{Predictions from gamlss2 Models}

\description{
Methods for \pkg{gamlss2} model objects that predict response distributions,
their parameters, additive predictors, individual terms, or response means.
}

\usage{
\method{predict}{gamlss2}(object, model = NULL, newdata = NULL,
  type = c("distribution", "parameter", "link", "response", "terms"), terms = NULL,
  se.fit = FALSE, drop = TRUE, ...,
  level = NULL, interval = c("none", "wald"), interval.cache = NULL)
}

\arguments{
  \item{object}{A model object of class \code{\link{gamlss2}}.}
  \item{model}{Character or numeric. Which distribution parameter(s) should be predicted?
    Names are family-specific, e.g., \code{"mu"}, \code{"sigma"}, or
    \code{"shape"}, \code{"rate"}. By default, predictions are returned for
    all model parts.}
  \item{newdata}{Data frame. Optionally, a new data frame in which to
    look for variables with which to predict. If omitted, the original
    observations are used.}
  \item{type}{Character. Which type of prediction should be computed? The
    default \code{"distribution"} returns \pkg{distributions3} objects. If
    \code{type} is omitted and \code{model}, \code{parameter}, or \code{what}
    is specified, \code{"parameter"} is used.
    Use \code{"link"} for additive predictors on the link scale,
    \code{"parameter"} for distribution parameters on their natural scale,
    \code{"terms"} for individual additive terms, and \code{"response"} for
    the mean of the fitted response distribution.}
  \item{terms}{Character. Which of the terms in the additive predictor(s)
    should be included? By default all terms are included.}
  \item{se.fit}{Logical. Should standard errors for the predictions be included?
    Standard errors are computed by simulating from the approximate multivariate normal
    distribution defined by the complete joint covariance of the maximum likelihood estimates.
    Local term-wise covariance is never used for these draws. The number of simulations is controlled
    by the argument \code{R}, which defaults to \code{R = 200}, and can be passed via \code{...}.}
  \item{drop}{Logical. Should predictions be simplified when possible?}
  \item{level}{\code{NULL}, or a numeric vector of confidence levels for
    Wald bands, with every value strictly between zero and one.}
  \item{interval}{Character. Use \code{"wald"} for fast confidence bands
    from the joint penalized information. The default \code{"none"} retains
    ordinary prediction behavior. Supplying \code{level} also selects Wald
    bands.}
  \item{interval.cache}{Optional joint covariance cache obtained from the
    \code{"interval.cache"} attribute of a previous Wald prediction for the
    same, unchanged model. Reusing it avoids rebuilding the joint covariance
    for subsequent predictions.}
  \item{...}{Further control arguments. Supported aliases are \code{what} and
    \code{parameter} for \code{model}. Additional arguments are used for
    simulation-based summaries, including \code{FUN}, \code{R}, \code{seed},
    \code{burnin}, and \code{nogrep}.}
}

\details{
Predictions can be computed for the original observations or for \code{newdata}.
By default, response distributions are returned. Depending on \code{type},
results can instead be returned as distribution parameters, on the link scale,
as individual model terms, or as the fitted response mean.

If \code{se.fit = TRUE}, standard errors are obtained from coefficient draws.
Existing draws stored in the fitted object are used when available; otherwise,
draws are generated from the complete joint covariance of the fitted
coefficients, including cross-term and cross-parameter dependence.

If \code{level} is specified, approximate pointwise confidence intervals and
standard errors are computed analytically, without simulation. Intervals on
transformed parameter or response scales are transformed as needed and may be
asymmetric. These are confidence intervals for the fitted model, not prediction
intervals for future observations, and they do not include smoothing-parameter
or model-selection uncertainty.

Wald intervals currently support maximum-likelihood fits with linear terms and
coefficient-linear \pkg{mgcv} smooths. They are not available for Bayesian fits
or nonlinear custom special terms. Wald intervals cannot be
combined with simulation controls such as \code{FUN}, \code{R}, \code{seed}, or
\code{burnin}. The optional \code{interval.cache} can be reused for repeated
predictions from the same unchanged model.

For \code{type = "terms"}, term names are partially matched by default. Use
\code{nogrep = TRUE} in \code{...} for exact matching.

See also \code{\link{prodist.gamlss2}} for explicitly obtaining
\pkg{distributions3} objects.
}

\value{
For \code{type = "distribution"}, a distribution object is returned. With
coefficient draws and \code{FUN = identity}, a list of distribution objects is
returned; the family is stored once for the complete list and added when an
element is extracted. For other types, \code{drop = TRUE} simplifies the result
when possible. Term predictions are data frames or, for multiple parameters,
lists of matrices. With \code{se.fit = TRUE}, fitted values and simulation-based
standard errors are returned.

When \code{level} is specified, a list with elements \code{fit}, \code{se.fit},
\code{lower}, and \code{upper} is returned. For one confidence level,
\code{lower} and \code{upper} follow the structure of the ordinary point
prediction, including \code{model}, \code{terms}, and \code{drop} selection.
For multiple confidence levels, each is a list named by confidence percentage
(for example, \code{"80\%"} and \code{"95\%"}), whose elements have that same
prediction structure. The \code{fit} element is the ordinary, non-simulated
prediction and \code{se.fit} is shared by all intervals. Attributes record the
confidence levels, method, and reusable \code{interval.cache}. The cache stores
the full joint covariance and fitted inferential state and can be reused only
with that unchanged model.

If the joint observed information is indefinite or a required coefficient is
not estimable, Wald prediction and coefficient simulation report the problem;
they do not silently substitute a working-information approximation.
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

## model selects parameter predictions if type is omitted
predict(m, newdata = nd, model = "sigma")
predict(m, newdata = nd, model = "sigma", drop = FALSE)

## individual terms in additive predictor(s)
predict(m, newdata = nd, type = "terms", model = "sigma")
predict(m, newdata = nd, type = "terms", model = "sigma", terms = "s(speed)")

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
q <- sapply(d, quantile, probs = 0.9)

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
