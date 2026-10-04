\name{JR}
\alias{JR}
\alias{jr}

\title{Joint REML Estimation}

\description{Fit a GAMLSS model by joint Laplace approximate REML.}

\usage{
## optimzer function
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
  If joint REML fails, the input fit is retained and \code{reml.failure}
  records the reason. Check \code{isTRUE(fit$jr$converged)} before using
  joint uncertainty estimates.

  The following controls can be supplied through
  \code{\link{gamlss2_control}} or as named arguments to \code{jr()}:
  \itemize{
    \item \code{jr.inner}: Function, \code{\link{RS}} or \code{\link{CG}},
      used to obtain the initial fit. Default is \code{RS}.
    \item \code{jr.outer.maxit}: Integer, the maximum number of outer
      smoothing parameter iterations. Default is \code{30}.
    \item \code{jr.outer.maxeval}: Integer, the maximum number of outer
      objective evaluations. Default is \code{150}.
  }
  The initial RS/CG fit uses \code{eps} and \code{maxit}, as documented in
  \code{\link{RS}}. The current direct coefficient solver does not use
  \code{jr.eps} or \code{jr.maxit}.

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
