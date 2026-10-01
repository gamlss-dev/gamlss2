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
  method = c("joint", "working", "numeric"), unconditional = FALSE,
  sandwich = FALSE, .details = FALSE, \dots)
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
  \item{sandwich}{Logical. If \code{TRUE}, use the post-fit estimating-equation
    sandwich covariance for an ordinary RS fit with local REML smooths.
    This requires \code{method = "joint"} and \code{unconditional = FALSE}.}
  \item{.details}{Internal diagnostic for the sandwich calculation. Ordinary
    calls should leave this at \code{FALSE}.}
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

With \code{sandwich = TRUE}, the coefficient score and the actual local REML
fixed-point equations form a joint estimating system for coefficients and log
smoothing parameters. Its bread uses the same joint penalized coefficient
information as the default method. Its meat is the outer product of
observation-wise derivatives with respect to case multipliers. Prior weights
are treated as fixed multipliers on independent rows; this is not a
frequency-replication interpretation. The quadratic penalty is deterministic
and is not included in the stochastic meat. The full joint coefficient and
smoothing-parameter covariance is assembled internally, while \code{vcov()}
returns only the coefficient block. No coefficient or smooth refit is performed.
This method requires a converged ordinary RS fit and local REML smoothing
estimates. An estimated smoothing parameter clipped at an RS limit is held
fixed when its local update before clipping points beyond that limit. A
single-penalty smooth whose local REML variance floor is active uses the
corresponding floored update equation. A user-fixed nonnegative smoothing
parameter is held fixed regardless of its
size, including values above the RS estimation cap. A stored estimated
parameter whose local update points inward is inconsistent with the final
fit. Indeterminate transition cases, non-REML smoothing selection, selected
ridge or adaptive penalties, and non-estimable coefficients are rejected.
It uses the final working equations and is a local fixed-point approximation
when RS has stopped at a finite convergence tolerance.

By default, the result is conditional on the fitted smoothing parameters.
For fits estimated with \code{optimizer = JR},
\code{unconditional = TRUE} uses the Hessian of the profiled joint REML
criterion stored by \code{\link{JR}}. The outer optimization includes
coefficient optimization at fixed smoothing parameters; no refits are required
when \code{vcov()} is called. A successful \code{JR} fit also stores the
conditional coefficient covariance from its direct likelihood mode. If the
joint REML optimum or its curvature cannot be
verified, unconditional covariance is unavailable. Curvature directions too
flat to distinguish from numerical error are excluded from the inverse for
variation of the coefficient mean. Their contribution to variation of the
coefficient covariance factor is bounded, as in \pkg{mgcv}.

For other optimizers, \code{unconditional = TRUE} uses a local quadratic
post-fit approximation to joint marginal smoothness curvature. The current
term-wise RS/CG smoothness selection does not provide the global REML Hessian
used by \pkg{mgcv}, so this post-fit correction is not a validated covariance
of the complete RS/CG estimator. Following \pkg{mgcv}, nonpositive directions
of the local curvature do not contribute to variation in the coefficient
mean. A small regularization bounds their contribution to variation in the
coefficient covariance factor. Adaptive smooth penalties, selected ridge
penalties, and non-estimable coefficients cause an explicit error. Fixed
smoothing parameters receive no correction. Model-selection uncertainty is
not included.

Special terms that have no verifiable coefficient-linear design are rejected
explicitly. Objects fitted with \code{light = TRUE} do not retain the final
predictors and must be refitted for joint covariance. Aliased or numerically
unidentified coefficient variances are returned as \code{NA}; no arbitrary
inferential ridge is added.

The ordinary \code{summary()} and effect-result paths use the covariance stored
with each fitted model term by default, which is substantially faster for large
models. Use \code{summary(object, method = "joint")} or
\code{results(object, method = "joint")} to request the complete joint
information. Use \code{summary(object, sandwich = TRUE)} to report standard
errors and coefficient tests from the RS sandwich covariance, or
\code{summary(object, unconditional = TRUE)} for the existing smoothness
correction. These options select the joint method unless another method is
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

## frequentist covariance for a regular RS/local-REML fit
Vs <- vcov(m, full = TRUE, sandwich = TRUE)
}}

\keyword{models}
\keyword{regression}
