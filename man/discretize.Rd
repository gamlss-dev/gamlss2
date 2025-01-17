\name{discretize}
\alias{discretize}

\title{Discretize Continuous Distribution Family for Count Regression Models}

\description{
  This function takes any continuous distribution family object and discretizes it,
  enabling it to be used for the estimation of count regression models. The discretized family
  can then be used in \code{\link{gamlss2}} models that deal with count data.
}
\usage{
discretize(family = NO)
}

\arguments{
  \item{family}{A continuous distribution family object. The family will be discretized
    for modeling count data, where the distribution is adapted for count outcomes.}
}

\value{
  Returns an object of class \code{"gamlss2.family"}, which is a discretized version of
  the input continuous family object, suitable for use in \code{\link{gamlss2}}
  models for count data.
}

\details{
  The function discretizes a continuous distribution family by converting its
  cumulative distribution function (CDF) into a probability mass function (PMF).
  This is done by computing the difference between the CDF evaluated at adjacent points.
  The resulting discretized distribution can be used in count regression models to estimate
  the relationship between count data and explanatory variables.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2.family}}
}

\examples{
## Simulate count data using the Poisson distribution.
set.seed(111)
y <- rpois(1000, lambda = 10)

## Create a discretized family using the BCT distribution (with log link for mu).
fam <- discretize(family = BCT(mu.link = "log"))

## Fit a count regression model using the discretized family.
fit_family(y, family = fam)
}

\keyword{utilities}
\keyword{count data}
\keyword{regression}
\keyword{discretization}

