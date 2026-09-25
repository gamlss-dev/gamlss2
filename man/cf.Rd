\name{cf}
\alias{cf}
\alias{special_fit.cf}
\alias{special_predict.cf.fitted}

\title{Conditional Inference Forest Model Terms}

\description{
  Constructs a conditional inference forest model term using \pkg{partykit}.
}

\usage{cf(formula, ...)}

\arguments{
  \item{formula}{A vector or formula specifiyng the covariates that should
    be used for modeling.}
  \item{\dots}{Control arguments passed to \code{\link[partykit]{cforest}}.}
}

\details{
  The forest is refitted with the current working response and weights
  during backfitting. The fitted contribution is centered for identifiability.
}

\value{
  An object of classes \code{"special"} and \code{"cf"}.
}

\seealso{
  \code{\link{special_terms}}
}

\examples{\dontrun{
data("mcycle", package = "MASS")

set.seed(6020)

m <- gamlss2(
  accel ~ cf(times) | .,
  data = mcycle
)

plot(m)

## predicted distributions and quantiles
predict(m, newdata = mcycle[1:5, , drop = FALSE])
pq <- quantile(m, probs = c(0.025, 0.5, 0.975))

plot(accel ~ times, data = mcycle,
  ylim = range(mcycle$accel, pq))
matlines(mcycle$times, pq, lty = 1, lwd = 2, col = 4)
}}

\keyword{models}
