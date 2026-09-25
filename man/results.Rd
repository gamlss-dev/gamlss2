\name{results}
\alias{results}
\alias{results.gamlss2}

\title{Extract Results from GAMLSS Models}

\description{
  Extracts partial effects and corresponding uncertainty intervals from a
  fitted \code{gamlss2} model. The results are mainly used for plotting.
}

\usage{
results(x, data = NULL, ...)

\method{results}{gamlss2}(x, data = NULL, ...)
}

\arguments{
  \item{x}{An object of class \code{"gamlss2"}.}
  \item{data}{An optional data frame used to construct the effects. By default,
    the model frame of \code{x} is used.}
  \item{\dots}{Further arguments controlling uncertainty intervals, including
    \code{interval}, \code{level}, \code{method}, and \code{nsim}.}
}

\value{
  A list containing the partial effects. Each effect is stored as a data frame
  with the corresponding covariate values and the columns \code{fit},
  \code{lower}, and \code{upper}.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{plot.gamlss2}}
}

\examples{\dontrun{## estimate model
m <- gamlss2(dist ~ s(speed), data = cars, family = GA)

## extract results
r <- results(m)
names(r$effects)
}}

\keyword{models}
\keyword{regression}
