\name{ct}
\alias{ct}
\alias{special_fit.ct}
\alias{special_predict.ct.fitted}

\title{Conditional Inference Tree Model Terms}

\description{
  Constructs a conditional inference tree model term using \pkg{partykit}.
}

\usage{
ct(formula, ...)
}
\details{The tree is refitted with the current working response and weights during backfitting. The fitted contribution is centered for identifiability.}
\value{An object of classes \code{"special"} and \code{"ct"}.}
\seealso{
  \code{\link{special_terms}}
}

\examples{\dontrun{
data("mcycle", package = "MASS")

m <- gamlss2(
  accel ~ ct(times,minsplit=5,alpha=0.9) | .,
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

