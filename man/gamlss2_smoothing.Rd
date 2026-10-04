\name{gamlss2_smoothing}
\alias{gamlss2_smoothing}

\title{Smoothing Parameter Selection}

\description{
  Smoothing parameter selection for smooth terms in models fitted with
  \code{\link{gamlss2}}.
}

\details{
  For \pkg{mgcv}-based smooth terms, smoothing parameters are selected by
  local REML updates by default (\code{criterion = "reml"}). These updates
  use the current working response and weights. With one penalty matrix this
  is the usual Schall update; with multiple penalty matrices all smoothing
  parameters are updated together using the extended Fellner--Schall trace
  equations. This is a local working-model REML update. Explicit choices such
  as \code{"gcv"}, \code{"aic"}, \code{"aicc"}, \code{"bic"}, and
  \code{"ncv"} retain their respective selection methods.

  These local updates are used by the \code{\link{RS}} and
  \code{\link{CG}} optimizers. With \code{optimizer = JR}, coefficients
  and smoothing parameters are estimated jointly by Laplace approximate
  REML; see \code{\link{JR}}.

  Neighbourhood cross-validation can be selected with \code{criterion = "ncv"}.
  If \code{ncv} is omitted, singleton neighbourhoods are used, corresponding to
  leave-one-observation-out cross-validation. For ordered data, a numeric
  \code{ncv = k} uses the observation together with up to \code{k} neighbours
  on either side. Larger neighbourhoods generally provide more smoothing and
  can account for short-range dependence, at increased computational cost.
  For ordered data, small values such as \code{ncv = 1} or \code{ncv = 2}
  are often preferable to the singleton default. Custom neighbourhoods can
  also be supplied as a list through a smooth term's \code{xt} argument, for
  example \code{s(x, xt = list(ncv = neighbourhoods))}, where
  \code{neighbourhoods} has one integer index vector per observation and each
  vector contains its target observation.

  All methods return coefficients and covariance at the reported smoothing parameters.

  For the \code{RS} and \code{CG} optimizers, argument
  \code{demmler.reinsch} controls the use of a weighted Demmler--Reinsch
  reparameterization when selecting the smoothing parameter of an eligible
  single-penalty smooth, including multivariate smooths. The default,
  \code{"auto"}, uses the native reparameterization for local REML updates,
  and for quadratic criteria after the initial outer iteration. Set
  \code{demmler.reinsch = TRUE} to force it for eligible updates, or
  \code{FALSE} to disable it. For quadratic criteria, the initial outer
  iteration continues to use the direct
  solver; with numerical gradients or \code{native.wfit = FALSE}, automatic
  reparameterization is delayed until the current search reaches the edge of
  its one-decade search window.

  Within an \code{RS} or \code{CG} fit, the reparameterization is cached
  while the design matrix, weights, penalty, and binning are unchanged. A new
  working response only requires updating its projection into the cached basis.
  The declared penalty rank keeps the null space unpenalized even at large
  smoothing parameters.

  Analytic gradients are used by default for quadratic smoothing criteria
  (GCV, AIC, GAIC, AICc, and BIC), including smooths with multiple penalties.
  Set \code{analytic.gradient = FALSE} to use numerical gradients. Criteria
  using \code{logLik = TRUE} continue to use numerical gradients, with native
  direct solves and family likelihood evaluation in R. Nearly exact fits use
  observation-space residuals to avoid cancellation in the quadratic RSS.

  Smooths with fixed smoothing parameters and multiple-penalty smooths use
  the direct solver. For local REML, the common range of the penalty matrices
  is cached and the coupled update is evaluated in native code. Single-penalty
  updates use the constrained penalty null space and the Demmler--Reinsch basis
  when available.
}

\references{
  Demmler, A. and Reinsch, C. (1975). Oscillation Matrices with Spline
  Smoothing. \emph{Numerische Mathematik}, 24(5), 375--382.
  \doi{10.1007/BF01437406}

  Wood, S. N. (2026). On Neighbourhood Cross Validation. \emph{Statistical
  Science}, in press. \doi{10.48550/arXiv.2404.16490}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2_control}},
  \code{\link{RS}}, \code{\link{CG}}, \code{\link{JR}},
  \code{\link{pb}}
}

\examples{\dontrun{## specify the model formula
f <- dist ~ s(speed) | s(speed)

## estimate models with different smoothing criteria
m0 <- gamlss2(f, data = cars, family = GA, criterion = "reml")
m1 <- gamlss2(f, data = cars, family = GA, criterion = "aic")
m2 <- gamlss2(f, data = cars, family = GA, criterion = "bic")
m3 <- gamlss2(f, data = cars, family = GA, criterion = "gaic")
m4 <- gamlss2(f, data = cars, family = GA, criterion = "gaic", K = 4)
m5 <- gamlss2(f, data = cars, family = GA, criterion = "ncv")
m6 <- gamlss2(f, data = cars, family = GA, criterion = "ncv", ncv = 2)

plot(c(m0, m1, m2, m3, m4, m5, m6))
}}

\keyword{models}
\keyword{regression}
