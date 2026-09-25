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

\examples{\dontrun{## air quality data
air <- subset(airquality, !is.na(Ozone))

set.seed(6020)
m <- gamlss2(
  Ozone ~ cf(Temp, minsplit = 2L, minbucket = 2L, alpha = 0.1) | .,
  data = air,
  family = GA
)

plot(m)

## predicted distributions and quantiles
d <- predict(m, newdata = air)
head(d, 3)
pq <- quantile(d, probs = c(0.025, 0.5, 0.975))

plot(Ozone ~ Temp, data = air,
  ylim = range(air$Ozone, pq))
i <- order(air$Temp)
matlines(air$Temp[i], pq[i, ], lty = 1, lwd = 2, col = 4)
}}

\keyword{models}

