\name{OL}
\alias{OL}
\alias{ologit}

\title{Ordered Logistic Family for Ordinal Regression}

\description{
  Defines the ordered logistic (cumulative logit) family for modeling ordinal response variables
  within the \code{gamlss2} framework. This implementation supports flexible modeling of the
  location and threshold (cutpoint) parameters, including effects of covariates.
}

\usage{
OL(k)

ologit(k)
}

\arguments{
  \item{k}{An integer specifying the number of response categories. Must be \code{k >= 2}.}
}

\details{
  This family implements a cumulative logit model for ordinal responses with \code{k}
  ordered categories. The linear predictor models a latent location parameter, and the
  cutpoints between response categories are parameterized via a monotonic transformation:
  \itemize{
    \item The first cutpoint is modeled directly (\code{theta1}).
    \item The remaining cutpoints are expressed as \code{theta1 + exp(delta_j)} for \code{j = 2, ..., k - 1}, 
      ensuring that the thresholds remain ordered.
  }

  The \code{OL()} family supports modeling the location and threshold differences (\code{delta_j}) 
  as functions of covariates using additive predictors in \code{\link{gamlss2}} via the \code{"|"} formula interface.

  The family returns an object of class \code{"gamlss2.family"}, which includes methods for evaluating 
  the log-likelihood, simulating from the model, and computing predicted probabilities.

  \code{ologit()} is retained as a compatibility alias for \code{OL()}.
}

\value{
  A \code{"gamlss2.family"} object to be used with \code{\link{gamlss2}}.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2.family}}, \code{\link[MASS]{polr}}
}

\examples{
## example using the housing data from the MASS package:
library("MASS")

## fit standard cumulative logit model using polr()
m1 <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
summary(m1)

## convert response to integer for use with gamlss2
housing$Satint <- as.integer(housing$Sat)

## fit equivalent model using gamlss2
m2 <- gamlss2(Satint ~ Infl + Type + Cont,
  data = housing, weights = Freq,
  family = OL(k = 3))
summary(m2)

## compare coefficients
coef(m1)
coef(m2)

## predict class probabilities
pm1 <- predict(m1, type = "p")
pm2 <- predict(m2)
pm2 <- family(m2)$probabilities(pm2)

print(head(pm1))
print(head(pm2))
}

\keyword{regression}
\keyword{ordinal}
\keyword{models}
\keyword{distribution}
