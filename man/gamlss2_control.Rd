\name{gamlss2_control} 
\alias{gamlss2_control}

\title{Control Parameters}

\description{
  Control parameters for fitting GAMLSS models with \code{\link{gamlss2}}.
}

\usage{
gamlss2_control(optimizer = RS, trace = TRUE,
  flush = TRUE, light = FALSE, expand = TRUE,
  model = TRUE, x = TRUE, y = TRUE,
  fixed = FALSE, \dots)
}

\arguments{
  \item{optimizer}{Function, the optimizer to be used for fitting.}
  \item{trace}{Logical, should information be printed while the algorithm is running?}
  \item{flush}{Logical, whether to use \code{\link[utils]{flush.console}} to display
    current output in the console.}
  \item{light}{Logical, if set to \code{light = TRUE}, no model frame, response,
    model matrix and other design matrices will be part of the return value.}
  \item{expand}{Logical, if fewer formulas are supplied than there are parameters
    of the distribution, should intercept-only formulas be added automatically?}
  \item{model}{Logical, should the \emph{model frame} be included as a component of the
    returned object.}
  \item{x}{Logical, indicating whether the model matrix should be included
    as a component of the returned object.}
  \item{y}{Logical, should the response be included as a component of the returned object.}
  \item{fixed}{Named logical vector indicating which parameters should be fixed during estimation.
    See \code{\link{gamlss2_start}} for examples.}
  \item{\dots}{Further control parameters to be included in the return value, for example
    parameters used by the optimizer function \code{\link{RS}}.}
}

\details{
  The set of control parameters can be extended. For example, if a different
  optimizer is used, newly specified control parameters are passed on to that
  optimizer automatically.

  For the \code{RS} and \code{CG} optimizers, argument
  \code{demmler.reinsch} controls the use of a weighted Demmler--Reinsch
  reparameterization when selecting the smoothing parameter of an eligible
  single-penalty smooth, including multivariate smooths. The default,
  \code{"auto"}, uses the native reparameterization for quadratic criteria
  after the initial outer iteration. Set \code{demmler.reinsch = TRUE} to
  force it for eligible updates, or \code{FALSE} to disable it. The initial
  outer iteration continues to use the direct solver. With numerical gradients
  or \code{native.wfit = FALSE}, automatic reparameterization is delayed until
  the current search reaches the edge of its one-decade search window.

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

  Smooths with fixed smoothing parameters or multiple penalties, including
  adaptive and tensor-product smooths, use the direct solver. For larger
  multiple-penalty smooths, a cached factor of the weighted crossproduct
  reduces the cost of EDF derivatives. Local ML updates use the constrained
  penalty null space and return coefficients and covariance at the reported
  smoothing parameter. Native single-penalty ML updates compute RSS and the
  penalty quadratic in the DR basis, with observation-space residuals for
  nearly exact fits.
} 

\value{
  A list with the arguments specified.
}

\seealso{
  \code{\link{RS}}, \code{\link{gamlss2}}, \code{\link{gamlss2_start}}
}

\examples{\dontrun{## specify the model formula
f <- dist ~ s(speed) | s(speed)

## estimate model with different step length
## control in the RS algorithm
m1 <- gamlss2(f, data = cars, family = GA, step = 1)
m2 <- gamlss2(f, data = cars, family = GA, step = 0.1)
}}

\keyword{models}
\keyword{regression}
