\name{MN}
\alias{MN}

\title{Multinomial Logit Family}

\description{
  The \code{MN()} family implements a multinomial logit model for
  categorical responses with \eqn{k \ge 2} unordered categories.
}

\usage{
MN(k)
}

\arguments{
  \item{k}{An integer specifying the number of response categories.}
}

\details{
  The response must be a factor with exactly \code{k} levels.
  The first level is treated as the reference category.

  For categories \eqn{j = 2, \ldots, k}, each predictor models the
  log-odds relative to the reference category.

  The family provides analytical first and second derivatives of the
  log-likelihood with respect to the predictors and a
  \code{probabilities()} method for fitted category probabilities.
}

\value{
  A \code{"gamlss2.family"} object to be used with \code{\link{gamlss2}}.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2.family}}
}

\examples{\dontrun{data("marital.nz", package = "VGAM")

## multinomial model
m <- gamlss2(mstatus ~ s(age) | . | ., data = marital.nz, family = MN(4))

## plot effects
plot(m)

## predict probabilities
par <- predict(m)
probs <- family(m)$probabilities(par)
names(probs) <- levels(marital.nz$mstatus)

## plot
i <- order(marital.nz$age)
matplot(marital.nz$age[i], probs[i, ], type = "l", lty = 1, lwd = 2,
  xlab = "Age", ylab = "Probabilities")
legend("center", names(probs), lty = 1, lwd = 2, col = 1:4)
}}

\keyword{regression}
\keyword{categorical}
\keyword{multinomial}
\keyword{distribution}
