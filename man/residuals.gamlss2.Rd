\name{residuals.gamlss2}
\alias{residuals.gamlss2}
\alias{plot.gamlss2.residuals}

\title{Residuals for GAMLSS Models}

\description{
  Extracts residuals from fitted \code{gamlss2} objects. By default,
  normalized quantile residuals are returned. Response residuals can be
  obtained as observed responses minus fitted distribution means. A plotting
  method for residual diagnostic plots is also provided.
}

\usage{
\method{residuals}{gamlss2}(object,
  type = c("quantile", "response", "parameter"), newdata = NULL, \dots)

\method{plot}{gamlss2.residuals}(x, which = NULL, spar = TRUE, \dots)
}

\arguments{
  \item{object}{An object of class \code{"gamlss2"}.}
  \item{type}{Character specifying the type of residuals to be computed.
    The default \code{"quantile"} computes normalized quantile residuals.
    \code{"response"} computes observed responses minus fitted distribution
    means. \code{"parameter"} is reserved for parameter-specific residuals
    and is currently not implemented.}
  \item{newdata}{An optional data frame in which to evaluate the residuals.
    If omitted, residuals are evaluated for the data used to fit the model.
    When supplied, \code{newdata} must also contain the response variable.}
  \item{\dots}{Additional arguments passed on to methods, family-specific
    residual functions, or \code{\link{predict.gamlss2}}.}
  \item{x}{An object of class \code{"gamlss2.residuals"}, typically produced
    by \code{residuals.gamlss2}.}
  \item{which}{Character or integer selecting the diagnostic plot(s).
    Available plots are \code{"hist-resid"} for a histogram with kernel
    density estimate, \code{"qq-resid"} for a normal quantile-quantile plot,
    and \code{"wp-resid"} for a worm plot. If omitted, all three plots are
    produced. Integers \code{1}, \code{2}, and \code{3} select these plots in
    the same order.}
  \item{spar}{Logical. Should graphical parameters be set automatically for
    multiple plots?}
}

\details{
  Normalized quantile residuals are obtained by evaluating the fitted
  distribution function at the observed responses and transforming the resulting
  probabilities to the standard normal scale. For continuous response
  distributions this is
  \deqn{r_i = \Phi^{-1}\{F(y_i | \hat\theta_i)\},}
  where \eqn{F()} is the fitted cumulative distribution function and
  \eqn{\hat\theta_i} denotes the fitted distribution parameters for observation
  \eqn{i}.

  For discrete response distributions, randomized quantile residuals are used:
  a uniform random probability is drawn between the fitted distribution
  probabilities just below and at the observed response before applying the
  standard normal quantile function. Thus, repeated calls may differ for
  discrete families unless the random-number seed is fixed.

  If the family object supplies a family-specific residual function, this is
  used. Otherwise, quantile residuals are computed from the family cumulative
  distribution function. If no cumulative distribution function is available,
  the function falls back to response residuals and issues a warning.

  Response residuals are computed as
  \deqn{y_i - \mu(\hat\theta_i),}
  where \eqn{\mu()} is the family mean function if available. Otherwise the
  first distribution parameter is used as the fitted mean.

  Non-finite residuals are replaced by \code{NA}. The plotting method provides
  simple residual diagnostics for checking whether quantile residuals are
  approximately standard normal.
}

\value{
  A numeric vector of residuals with class \code{"gamlss2.residuals"} and an
  attribute \code{"type"} indicating the residual type. The plotting method is
  called for its side effect and returns \code{NULL} invisibly.
}

\references{
  Dunn PK, Smyth GK (1996). \dQuote{Randomized Quantile Residuals.}
  \emph{Journal of Computational and Graphical Statistics}, \bold{5}(3),
  236--244. \doi{10.1080/10618600.1996.10474708}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{predict.gamlss2}},
  \code{\link{plot.gamlss2}}
}

\examples{
## data
data("cars", package = "datasets")

## fit heteroscedastic normal GAMLSS model
m <- gamlss2(dist ~ s(speed) | s(speed), data = cars, family = NO)

## plot estimated quantiles
p <- quantile(m, probs = c(0.05, 0.5, 0.95))
plot(dist ~ speed, data = cars, ylim = range(cars$dist, p))
matlines(cars$speed, p, lty = 1, col = 4, lwd = 2)

## normalized quantile residuals
rq <- residuals(m)
summary(rq)
attr(rq, "type")

## response residuals
rr <- residuals(m, type = "response")
summary(rr)

## residuals for new data
nd <- cars[1:10, ]
residuals(m, newdata = nd)

## diagnostic plots
plot(rq)
plot(rq, which = "qq-resid")
}

\keyword{models}
\keyword{regression}

