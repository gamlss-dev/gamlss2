\name{OP}
\alias{OP}
\alias{oprobit}

\title{Ordered Probit Family for Ordinal Regression}

\description{
  Defines the ordered probit family for ordinal responses in \code{gamlss2}.
}

\usage{
OP(k)

oprobit(k)
}

\arguments{
  \item{k}{A single integer specifying the number of response categories.
    Must be \code{k >= 2}.}
}

\details{
  The response must be coded as the integers \code{1, ..., k}, or as a
  factor with exactly these levels.

  The cumulative probabilities are
  \deqn{\Pr(Y > j) = \Phi(\mathrm{location} - \gamma_j), \quad j = 1, \ldots, k - 1,}
  where \eqn{\Phi} is the standard normal distribution function. As in
  \code{\link{OL}}, the cutpoints are parameterized by
  \eqn{\gamma_1 = \theta_1} and
  \eqn{\gamma_j = \gamma_{j-1} + \exp(\delta_j)} for
  \eqn{j = 2, \ldots, k - 1}, ensuring strictly ordered cutpoints.

  Each parameter can use an additive predictor through the \code{"|"}
  formula interface. The family provides category and transition
  probabilities, moments of numeric category labels, and randomized
  quantile residuals. \code{oprobit()} is an alias for \code{OP()}.

  Analytic scores, negative Hessians (including cross derivatives), and
  a joint working response and weight update are provided for estimation.
}

\value{
  A \code{"gamlss2.family"} object with parameters \code{location},
  \code{theta1}, and, for \code{k > 2}, \code{delta}\eqn{j} for
  \eqn{j = 2, \ldots, k - 1}, all with identity links.
}

\seealso{
  \code{\link{OL}}, \code{\link{gamlss2}}, \code{\link{gamlss2.family}},
  \code{\link[MASS]{polr}}
}

\examples{
## example using the housing data from the MASS package
data("housing", package = "MASS")
housing$Satint <- as.integer(housing$Sat)

## fit cumulative probit model
m <- gamlss2(Satint ~ Infl + Type + Cont,
  data = housing, weights = Freq, family = OP(k = 3))

## predict class probabilities
par <- predict(m, type = "parameter")
head(family(m)$probabilities(par))
}

\keyword{regression}
\keyword{ordinal}
\keyword{models}
\keyword{distribution}
