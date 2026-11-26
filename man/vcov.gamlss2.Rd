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
  \item{unconditional}{Logical. If \code{TRUE}, a correction for
    uncertainty in selected smoothness parameters is added. For \code{JR} fits this
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

Joint information is computed by differentiating the linked family scores.
Mixed derivatives are checked in both orientations. If finite differences
disagree, additional step sizes and Richardson extrapolation are used to
check numerical convergence; persistent inconsistencies raise an error.

For monotonic \code{\link{ms}} terms, inversion is restricted to coefficient
directions that preserve the fitted active constraints and term centering.
This local approximation conditions on the active constraint set and does not
include uncertainty about which constraints are active. Unconditional and
sandwich corrections are currently unavailable for these terms.

For \code{\link{fk}} terms, the identifying coefficient fixed to zero during
fitting is kept fixed in the covariance calculation. Inference conditions on
the selected knot locations and number of knots. Smoothness and sandwich
corrections are currently unavailable for these constrained representations.

The default covariance conditions on fitted smoothing parameters. Use
\code{unconditional = TRUE} to include their approximate estimation uncertainty.
Fixed smoothing parameters receive no correction, and model-selection
uncertainty is excluded. The smoothness correction does not support adaptive
smooth penalties, selected ridge penalties, or non-estimable coefficients.

Special terms that have no verifiable coefficient-linear design are rejected
explicitly. A stored term covariance alone is insufficient for adjusting
linear-coefficient covariance: the coefficient design, fitted quadratic
penalty (possibly zero), and identifying constraints must also be available.
Standard \pkg{mgcv} terms provide their design and penalty matrices; other
specials must supply an equivalent verifiable representation.
Objects fitted with \code{light = TRUE} do not retain the final
predictors and must be refitted for joint covariance. Aliased or numerically
unidentified coefficient variances are returned as \code{NA}.

The default \code{summary(object, method = "local")} and
\code{confint(object, method = "local")} extract the linear covariance
stored with the fitted model terms. For RS and CG fits this is computed once
after fitting using predictor-local working information, equivalent to
\code{vcov(object, method = "working")} for the linear coefficients.
All linear and smooth coefficients within each
distributional predictor enter the penalized inverse before extracting its
linear block. Thus random effects and smooths are not treated as known
offsets. This approximation omits cross-parameter curvature and conditions
on fitted smoothing parameters; it need not agree with the observed joint
covariance. Direct JR fits already supply marginal linear blocks from their
joint covariance, which are retained. Summaries and confidence intervals do
not reconstruct or refactor these systems. Light fits retain the same
corrected linear blocks. Models fitted before this correction must be refitted
to obtain corrected local inference in predictors containing smooth terms.
If a predictor contains unsupported specials, its local linear covariance is
unavailable and its standard errors are reported as \code{NA} with a warning.

Effect-result paths use the covariance stored with each fitted model term by
default. Use \code{summary(object, method = "joint")} or
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
