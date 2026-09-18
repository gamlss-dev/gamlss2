\name{tree}
\alias{tree}
\alias{special_fit.tree}
\alias{special_predict.tree.fitted}

\title{Regression Tree Model Terms}

\description{
  Constructs an \code{rpart} regression-tree model term for use with
  \code{gamlss2}.
}

\arguments{
  \item{formula}{A one-sided formula specifying the predictors in the tree.}
  \item{...}{Arguments passed to \code{\link[rpart]{rpart.control}}.}
}

\details{
  The term is refitted with \code{\link[rpart]{rpart}} to the current
  working response and weights at each backfitting update. Its contribution
  is centered before being returned to the additive predictor.
}

\value{An object of classes \code{"special"} and \code{"tree"}.}

\seealso{
  \code{\link{special_terms}}
}

\examples{\dontrun{
data("mcycle", package = "MASS")

m <- gamlss2(
  accel ~ tree(times,minsplit=5) | .,
  data = mcycle
)

plot(m)

predict(m, newdata = mcycle[1:5, , drop = FALSE])
pq <- quantile(m, probs = c(0.025, 0.5, 0.975))

plot(accel ~ times, data = mcycle,
  ylim = range(mcycle$accel, pq))
matlines(mcycle$times, pq, lty = 1, lwd = 2, col = 4)
}}

\keyword{models}

