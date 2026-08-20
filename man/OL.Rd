\name{OL}
\alias{OL}
\alias{ologit}

\title{Ordered Logistic Family for Ordinal Regression}

\description{
  Defines the ordered logistic (cumulative logit) family for modeling ordinal
  response variables within the \code{gamlss2} framework. The location and
  ordered cutpoints can be modeled as functions of covariates.
}

\usage{
OL(k)

ologit(k)
}

\arguments{
  \item{k}{A single integer specifying the number of response categories.
    Must be \code{k >= 2}.}
}

\details{
  This family implements a cumulative logit model for an ordinal response with
  \code{k} categories. The response must be coded as the integers
  \code{1, ..., k}, or as a factor with exactly these levels.

  The linear predictor models a latent \code{location} parameter. Let
  \eqn{\gamma_j} denote the cutpoint between categories \eqn{j} and
  \eqn{j + 1}. The cutpoints are parameterized as
  \deqn{\gamma_1 = \theta_1}
  and, recursively,
  \deqn{\gamma_j = \gamma_{j-1} + \exp(\delta_j),
    \quad j = 2, \ldots, k - 1.}
  Equivalently,
  \deqn{\gamma_j = \gamma_1 + \sum_{\ell=2}^j \exp(\delta_\ell),}
  which ensures strictly ordered cutpoints. In the family object,
  \eqn{\theta_1} is named \code{theta1} and \eqn{\delta_j} is named
  \code{delta}\eqn{j}. For \code{k = 2}, only \code{location} and
  \code{theta1} are required; there are no \code{delta} parameters.

  The cumulative probabilities are
  \deqn{\Pr(Y > j) = \mathrm{logit}^{-1}(\mathrm{location} - \gamma_j),
    \quad j = 1, \ldots, k - 1.}

  The \code{OL()} family supports modeling the location, first cutpoint, and
  cutpoint differences as functions of covariates using additive predictors
  in \code{\link{gamlss2}} via the \code{"|"} formula interface.

  The returned family object provides functions for evaluating probability
  masses and log-likelihoods, cumulative probabilities and quantiles,
  category and transition probabilities, and the mean and variance of the
  numeric category labels.

  \code{ologit()} is retained as a compatibility alias for \code{OL()}.
}

\value{
  A \code{"gamlss2.family"} object to be used with \code{\link{gamlss2}}.
  Its distribution parameters are \code{location}, \code{theta1}, and, when
  \code{k > 2}, one cutpoint-difference parameter \code{delta}\eqn{j} for
  each \eqn{j = 2, \ldots, k - 1}.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2.family}}, \code{\link[MASS]{polr}}
}

\examples{
## example using the housing data from the MASS package
library("MASS")

## fit standard cumulative logit model using polr()
m1 <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
summary(m1)

## OL() uses integer category labels 1, ..., k
housing$Satint <- as.integer(housing$Sat)

## fit the corresponding model using gamlss2
m2 <- gamlss2(Satint ~ Infl + Type + Cont,
  data = housing, weights = Freq,
  family = OL(k = 3))
summary(m2)

## compare coefficients
coef(m1)
coef(m2)

## predict class probabilities
pm1 <- predict(m1, type = "p")
pm2 <- family(m2)$probabilities(predict(m2))

print(head(pm1))
print(head(pm2))
}

\keyword{regression}
\keyword{ordinal}
\keyword{models}
\keyword{distribution}
