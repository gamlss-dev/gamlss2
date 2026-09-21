\name{vcov.gamlss2}
\alias{vcov.gamlss2}
\title{Joint Covariance Matrix for gamlss2 Models}

\description{
  Computes coefficient covariance from the joint penalized information of a
  fitted \code{gamlss2} model.
}

\usage{
\method{vcov}{gamlss2}(object,
  type = c("vcov", "cor", "se", "coef"), full = FALSE,
  method = c("joint", "working", "numeric"), \dots)
}

\arguments{
  \item{object}{A fitted \code{\link{gamlss2}} model.}
  \item{type}{The covariance matrix, correlation matrix, standard errors, or
    corresponding coefficient vector to return.}
  \item{full}{Logical. If \code{FALSE}, return ordinary linear coefficients.
    Their covariance is still extracted from the complete joint inverse. If
    \code{TRUE}, include supported smooth coefficients.}
  \item{method}{Information calculation. \code{"joint"} uses the raw joint
    observed information, \code{"working"} uses the final family working
    curvature with zero cross-parameter blocks, and \code{"numeric"} uses a
    slow coefficient-dimensional numerical Hessian for validation.}
  \item{\dots}{Control arguments passed to the numerical Hessian when
    \code{method = "numeric"}.}
}

\details{
The default covariance is the inverse of \eqn{H + P}, where \eqn{H} is the
joint likelihood information for every supported distributional predictor and
\eqn{P} is the complete penalty at the fitted smoothing and regularization
parameters. Thus it retains covariance among linear and smooth terms and among
distributional parameters. Prior observation weights and offsets in the final
linear predictors are included.

The result is conditional on the fitted smoothing parameters. It does not
include smoothing-parameter or model-selection uncertainty. Information and
its factorization are built lazily and cached for an unchanged fitted state.

Nonlinear special terms that have no coefficient-linear design are rejected
explicitly. Objects fitted with \code{light = TRUE} do not retain the final
predictors and must be refitted for joint covariance. Aliased or numerically
unidentified coefficient variances are returned as \code{NA}; no arbitrary
inferential ridge is added.

Effect plots and Wald predictions require a positive precision factor. If the
observed joint information is indefinite, Wald prediction warns and falls back
to \code{method = "working"}. Effect plots silently use the covariance stored
with each fitted term. Calling \code{vcov()} itself does not silently apply
either approximation.
}

\value{
A matrix for \code{type = "vcov"} or \code{"cor"}; a named numeric vector
for \code{type = "se"} or \code{"coef"}.
}

\seealso{
\code{\link{gamlss2}}, \code{\link{predict.gamlss2}}
}

\examples{\dontrun{
data("cars", package = "datasets")
m <- gamlss2(dist ~ s(speed) | s(speed), data = cars, family = NO)

## linear-coefficient covariance,
## marginalized over smooth coefficients.
vcov(m)

## complete linear and smooth covariance
V <- vcov(m, full = TRUE)
head(diag(V))
}}

\keyword{models}
\keyword{regression}
