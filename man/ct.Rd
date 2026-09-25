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

\arguments{
  \item{formula}{A vector or formula specifiyng the covariates that should
    be used for modeling.}
  \item{\dots}{Control arguments passed to \code{\link[partykit]{ctree}}.}
}

\details{
  The tree is refitted with the current working response and weights during 
  backfitting. The fitted contribution is centered for identifiability.
}

\value{
  An object of classes \code{"special"} and \code{"ct"}.
}

\seealso{
  \code{\link{special_terms}}
}

\examples{\dontrun{
m <- gamlss2(
  mag ~ ct(depth) | .,
  data = quakes,
  family = truncate_family(NO, lower = 3.99)
)

plot(m)

## predicted distributions and quantiles
d <- predict(m, newdata = quakes)
pq <- quantile(d, probs = c(0.025, 0.5, 0.975))

plot(mag ~ depth, data = quakes,
  ylim = range(quakes$mag, pq),
  col = adjustcolor(1, 0.3))
i <- order(quakes$depth)
matlines(quakes$depth[i], pq[i, ], lty = 1, lwd = 2, col = 4)
}}

\keyword{models}

