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
  method = c("joint", "working", "numeric"), unconditional = FALSE, \dots)
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
  \item{unconditional}{Logical. If \code{TRUE}, add a correction for
    uncertainty in selected smoothness parameters. For \code{JR} fits this
    uses the fitted joint REML curvature. Other fits use a post-fit local
    approximation. This requires \code{method = "joint"}.}
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

The default covariance conditions on fitted smoothing parameters. Use
\code{unconditional = TRUE} to include their approximate estimation uncertainty.
Fixed smoothing parameters receive no correction, and model-selection
uncertainty is excluded. Adaptive smooth penalties, selected ridge penalties,
and non-estimable coefficients are unsupported.

For \code{\link{JR}} fits, the correction uses the stored Hessian of the
profiled joint REML criterion and conditional coefficient covariance from the
direct likelihood mode; no refits are required. It is available only when the
joint REML optimum and its curvature can be verified. For other optimizers,
the correction uses a local quadratic post-fit approximation, rather than a
global REML Hessian. It is therefore not a validated covariance of the complete
RS/CG estimator.

Following \pkg{mgcv}, numerically flat directions for \code{JR}, or nonpositive
local curvature directions for other optimizers, are excluded from variation
in the coefficient mean. Their contribution to variation in the coefficient
covariance factor is bounded.

Special terms that have no verifiable coefficient-linear design are rejected
explicitly. Objects fitted with \code{light = TRUE} do not retain the final
predictors and must be refitted for joint covariance. Aliased or numerically
unidentified coefficient variances are returned as \code{NA}; no arbitrary
inferential ridge is added.

The ordinary \code{summary()} and effect-result paths use the covariance stored
with each fitted model term by default, which is substantially faster for large
models. Use \code{summary(object, method = "joint")} or
\code{results(object, method = "joint")} to request the complete joint
information. Use \code{summary(object, unconditional = TRUE)} for the smoothness
correction. This option selects the joint method unless another method is
specified explicitly. Simulation-based \code{predict()} calls always use the
complete joint covariance. No consumer silently substitutes \code{method = "working"}
when observed joint information is indefinite.
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

## include estimated smoothness uncertainty
Vc <- vcov(m, full = TRUE, unconditional = TRUE)
}}

\keyword{models}
\keyword{regression}
