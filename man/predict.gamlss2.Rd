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
  se.fit = FALSE, drop = TRUE, ...,
  level = NULL, interval.cache = NULL)
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
  \item{level}{\code{NULL}, or a numeric vector of confidence levels for
    Wald bands, with every value strictly between zero and one.}
  \item{interval.cache}{Optional information factor obtained from the
    \code{"interval.cache"} attribute of a previous Wald prediction for the
    same, unchanged model. Reusing it avoids rebuilding the joint information
    matrix for subsequent predictions.}
  \item{...}{Further control arguments. Supported aliases are \code{what} and
    \code{parameter} for \code{model}. Additional arguments are used for
    simulation-based summaries, including \code{FUN}, \code{R}, \code{seed},
    \code{burnin}, and \code{nogrep}.}
}

\details{
Predictions can be computed for the original observations or for \code{newdata}.
Depending on \code{type}, results are returned on the link scale, as distribution
parameters, as individual model terms, or as the fitted response mean.

If \code{se.fit = TRUE}, standard errors are obtained from coefficient draws.
Existing draws stored in the fitted object are used when available; otherwise,
draws are generated from the approximate distribution of the fitted coefficients.

If \code{level} is specified, approximate pointwise confidence intervals and
standard errors are computed analytically, without simulation. Intervals on
transformed parameter or response scales are transformed as needed and may be
asymmetric. These are confidence intervals for the fitted model, not prediction
intervals for future observations, and they do not include smoothing-parameter
or model-selection uncertainty.

Wald intervals currently support maximum-likelihood fits with linear terms and
standard \pkg{mgcv} smooths. They are not available for Bayesian fits, custom
special terms, or automatic term-selection penalties. Wald intervals cannot be
combined with simulation controls such as \code{FUN}, \code{R}, \code{seed}, or
\code{burnin}. The optional \code{interval.cache} can be reused for repeated
predictions from the same unchanged model.

For \code{type = "terms"}, term names are partially matched by default. Use
\code{nogrep = TRUE} in \code{...} for exact matching.

See also \code{\link{prodist.gamlss2}} for creating a full
\pkg{distributions3} object for moments, probabilities, quantiles, and random
numbers.
}

\value{
If \code{drop = FALSE}, a data frame is returned. If \code{drop = TRUE}
(the default), the result may be simplified to a numeric vector when possible.
For \code{type = "terms"}, a data frame is returned with one column per selected term.
Multi-parameter models return a list of term matrices.
For \code{se.fit = TRUE}, the returned object contains fitted values and
corresponding simulation-based standard errors.

When \code{level} is specified, a list with elements \code{fit}, \code{se.fit},
\code{lower}, and \code{upper} is returned. For one confidence level,
\code{lower} and \code{upper} follow the structure of the ordinary point
prediction, including \code{model}, \code{terms}, and \code{drop} selection.
For multiple confidence levels, each is a list named by confidence percentage
(for example, \code{"80\%"} and \code{"95\%"}), whose elements have that same
prediction structure. The \code{fit} element is the ordinary, non-simulated
prediction and \code{se.fit} is shared by all intervals. Attributes record the
confidence levels, method, and reusable \code{interval.cache}. The cache retains
a reference to the fitted model; reuse it only with that unchanged model.
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

## fast pointwise confidence bands, without sampling
bands <- predict(m, newdata = nd, level = c(0.8, 0.95))
bands$lower[["80\%"]]
bands$upper[["95\%"]]
}

\keyword{regression}
