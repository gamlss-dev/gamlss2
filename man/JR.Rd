\name{JR}
\alias{JR}
\alias{jr}

\title{Joint REML Estimation}

\description{Fit a GAMLSS model by joint Laplace approximate REML.}

\usage{
## optimizer function
JR(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)

## refit function
jr(object, trace = TRUE, ...)
}

\arguments{
  \item{x, y, specials, family, offsets, weights, start, xterms, sterms, control}{
    Model components passed by \code{gamlss2()}.}
  \item{object}{A fitted \code{gamlss2} model retaining its design matrices.}
  \item{trace}{Print the final joint REML result.}
  \item{\dots}{Named joint REML control overrides.}
}

\details{
  Use \code{optimizer = JR} to fit a model or \code{jr(fit)} to continue from
  an existing fit. The direct solver jointly estimates coefficients and
  smoothing parameters from cached designs and quadratic penalties. It
  requires a family density, a predictor-to-parameter map, and an
  identifiable finite joint mode.
  Selected penalties and fixed distributional parameters are unsupported.
  Duplicate unpenalized directions in overlapping smooths are removed from
  the joint solve. The stored coefficients and covariance use the
  minimum-norm representation in the original coefficient space, and
  effective degrees of freedom exclude these redundant directions.
  If joint REML fails, the input fit is retained and \code{reml.failure}
  records the reason. Check \code{isTRUE(fit$jr$converged)} before using
  joint uncertainty estimates.
  Power-exponential families may have a non-smooth REML profile when a
  fitted location coincides with an observation. Such fits retain RS/CG
  and report that smoothing-parameter curvature could not be verified.
  Fully penalized smooths that start nearly suppressed are also checked at
  less penalized starting values, to avoid stopping in a flat profile tail.
  If the outer optimizer reports false convergence, a stricter independent
  profile Newton check may verify the fit. In that case
  \code{fit$jr$outer.verified} is true; the original optimizer status is
  retained in \code{fit$jr$outer}. Iteration and evaluation limits still
  retain the initial fit.

  The following controls can be supplied through
  \code{\link{gamlss2_control}} or as named arguments to \code{jr()}:
  \itemize{
    \item \code{jr.inner}: Function, \code{\link{RS}} or \code{\link{CG}},
      used to obtain the initial fit. Default is \code{RS}.
    \item \code{jr.outer.maxit}: Integer, the maximum number of outer
      smoothing parameter iterations. Default is \code{60}.
    \item \code{jr.outer.maxeval}: Integer, the maximum number of outer
      objective evaluations. Default is \code{150}.
    \item \code{jr.eps}: Positive number, the relative coefficient step
      tolerance in the joint mode solve. Default is \code{1e-8}. The solver
      may refine further to stabilize power-exponential curvature.
    \item \code{jr.maxit}: Positive integer, the maximum number of
      coefficient iterations per objective evaluation. Default is \code{100}.
  }
  The initial RS/CG fit uses \code{eps} and \code{maxit}, as documented in
  \code{\link{RS}}. Fully and partially fixed smoothing parameters and
  unpenalized \code{fx = TRUE} smooths are supported. At extreme smoothing
  parameters, penalty roots and an equilibrated solve preserve accuracy
  in the unpenalized coefficient directions.

  For local REML updates with RS/CG and other smoothing parameter selection
  methods, see \code{\link{gamlss2_smoothing}}.
}

\value{
  A fitted \code{gamlss2} model with joint REML results in \code{fit$jr}.
}

\seealso{
  \code{\link{RS}}, \code{\link{gamlss2_control}},
  \code{\link{gamlss2_smoothing}}, \code{\link{vcov.gamlss2}}
}

\examples{\dontrun{
m <- gamlss2(dist ~ s(speed), data = cars, family = NO)
m <- jr(m)
vcov(m, unconditional = TRUE)
}}

\keyword{models}
