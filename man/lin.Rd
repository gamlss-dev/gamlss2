\name{lin}
\alias{lin}
\alias{ridge}
\alias{smooth.construct.lin.smooth.spec}
\alias{Predict.matrix.lin.effect}

\title{Bundled Linear and Ridge-Penalized Model Terms}

\description{
  Constructs a model-matrix block for use as a single model term in a
  \code{gamlss2} formula. Function \code{lin()} fits the block with only a
  negligible numerical penalty, whereas \code{ridge()} standardizes its
  columns and estimates a ridge penalty for the complete block.
}

\usage{
lin(x, ..., ridge = FALSE, scale = FALSE)

ridge(...)
}

\arguments{
  \item{x}{
    A one-sided formula or the first unquoted predictor expression. A formula
    is recommended when transformations, factors, or interactions are used,
    for example \code{~ x1 + x2}, \code{~ x1:x2}, or \code{~ 0 + f}.
  }
  \item{\dots}{
    In \code{lin()}, additional unquoted predictor expressions that are
    combined additively with \code{x} when \code{x} is not a formula. In
    \code{ridge()}, all arguments are forwarded to \code{lin()}; consequently,
    calls such as \code{ridge(x1, x2)} and \code{ridge(~ x1 + x2)} are both
    supported. For precise control of formula expansion, a one-sided formula
    is recommended.
  }
  \item{ridge}{
    Logical. Should the identity-penalty smoothing parameter be estimated?
    This is primarily a low-level argument; \code{ridge()} is the convenient
    interface for \code{ridge = TRUE} with scaling enabled.
  }
  \item{scale}{
    Logical. Should the columns of the constructed design matrix be centered
    and divided by their sample standard deviations? Columns with zero or
    non-finite standard deviation are left unchanged. Function \code{ridge()}
    always uses \code{scale = TRUE}.
  }
}

\details{
  \strong{Model-matrix construction.}
  The supplied predictors are converted to a design matrix using
  \code{\link[stats]{model.matrix}}. Any column named \code{"(Intercept)"}
  is removed because the distributional predictor has its own intercept.
  Apart from that removal, standard R formula rules and contrasts are used.
  Thus, \code{~ x1 * x2} expands to two main effects and their interaction,
  while \code{~ x1:x2} produces a pure interaction without adding the
  marginal terms. Pure interactions are therefore supported. For a factor
  \code{f}, \code{~ f} uses the current contrast coding (normally treatment
  coding), whereas \code{~ 0 + f} requests one indicator column per level.

  Supplying several predictors without a formula is a convenience shortcut:
  \code{lin(x1, x2)} is equivalent to \code{lin(~ x1 + x2)}, and similarly
  for \code{ridge()}. A formula should be used whenever operators such as
  \code{:}, \code{*}, \code{-}, or \code{0 +} need to be interpreted
  explicitly.

  \strong{Unpenalized linear blocks.}
  Function \code{lin()} fixes the identity-penalty multiplier at a negligible
  value. It therefore behaves like weighted linear regression for the current
  working response, while entering all columns as one model term rather than
  as separate ordinary linear terms. This can be useful for keeping a block of
  columns together and for using the same model-matrix and prediction machinery
  as penalized terms. Ordinary linear effects can still be specified directly
  in the surrounding \code{gamlss2} formula.

  \strong{Ridge-penalized blocks.}
  For \code{ridge()}, coefficients are estimated from the current weighted
  working response using an identity penalty. Ignoring the small numerical
  stabilization used by the fitting algorithm, the update minimizes

  \deqn{(z - X\beta)^\mathsf{T} W (z - X\beta) + \lambda\,\beta^\mathsf{T}\beta,}

  where \eqn{z} and \eqn{W} are the current working response and weights.
  A large \eqn{\lambda} shrinks the whole coefficient block toward zero.

  Before fitting, every nonconstant column of \eqn{X} is standardized using
  its training-data mean and standard deviation. The stored transformation is
  reused for prediction. Consequently, the identity penalty is applied on a
  comparable scale across numeric predictors, indicator columns, and
  interaction columns.

  The smoothing parameter \eqn{\lambda} is selected during backfitting using
  the smooth-term criterion supplied to \code{gamlss2()}. The default is
  corrected AIC (\code{criterion = "aicc"}); supported alternatives include
  \code{"gcv"}, \code{"aic"}, \code{"gaic"}, and \code{"bic"}. For
  example, use \code{gamlss2(..., criterion = "bic")} to select the ridge
  penalty by BIC. The previous estimate is used to warm-start subsequent
  backfitting updates.

  The effective degrees of freedom are computed from the penalized weighted
  hat matrix.
}

\value{
  \code{lin()} and \code{ridge()} return a smooth specification object of
  class \code{"lin.smooth.spec"}, which is completed to an object of class
  \code{"lin.effect"} during model setup.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{specials}},
  \code{\link[stats]{model.matrix}}, \code{\link[mgcv]{smooth.construct}}
}

\examples{\dontrun{
## longley data example with
## strongly correlated predictors
data("longley")

## fit all predictors linearly
m1 <- gamlss2(
  Employed ~ lin(~ GNP.deflator + GNP + Unemployed +
    Armed.Forces + Population + Year),
  data = longley,
  family = NO
)

## fit the same model using a ridge-penalized block,
## ridge() standardizes the model-matrix columns before applying the
## common ridge penalty.
m2 <- gamlss2(
  Employed ~ ridge(~ GNP.deflator + GNP + Unemployed +
    Armed.Forces + Population + Year),
  data = longley,
  family = NO
)

## Compare model summaries.
summary(m1)
summary(m2)

## extract coefficients
cm <- cbind(
  "linear" = specials(m1, model = "mu", elements = "coefficients"),
  "ridge" = specials(m2, model = "mu", elements = "coefficients")
)
print(cm)

## ridge terms combined with ordinary linear effects,
## Year is left unpenalized, while the remaining correlated
## economic variables are shrunk jointly
m3 <- gamlss2(
  Employed ~ Year +
    ridge(~ GNP.deflator + GNP + Unemployed +
      Armed.Forces + Population),
  data = longley,
  family = NO
)

summary(m3)

## interactions
m4 <- gamlss2(
  Employed ~ ridge(~ GNP + Unemployed + GNP:Unemployed),
  data = longley,
  family = NO
)

specials(m4, model = "mu", elements = "coefficients")
}}

\keyword{models}
\keyword{regression}
\keyword{smooth}

