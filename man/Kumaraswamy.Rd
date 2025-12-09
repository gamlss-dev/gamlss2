\name{Kumaraswamy}
\alias{Kumaraswamy}
\alias{KS}
\alias{shiftlog}

\title{Kumaraswamy Distribution}

\description{
  This function implements the two-parameter Kumaraswamy family for responses in (0, 1).
}

\usage{
## The Kumaraswamy family.
Kumaraswamy(a.link = shiftlog, b.link = shiftlog, ...)
KS(a.link = shiftlog, b.link = shiftlog, ...)

## The exp(x) + shift link specification.
shiftlog(shift = 1)
}

\arguments{
  \item{a.link}{Character or function, the link function to be used for parameter \code{a}.}
  \item{b.link}{Character or function, the link function to be used for parameter \code{b}.}
  \item{shift}{Numeric, the shift parameter to be used for the link.}
  \item{\dots}{Not used.}
}

\details{
The Kumaraswamy distribution is a continuous distribution defined on the interval (0, 1). The probability density function is

\deqn{\displaystyle f(y; a, b) = aby^{a - 1}(1 - y^a)^{b - 1}}

\eqn{y \in (0, 1)} is the response, \eqn{a} and \eqn{b} are non-negative parameters.

The \code{shiftlog} link function is given by:

\deqn{\displaystyle \exp(x) + 1}

This is the default, since the mode of the distribution is only defined for \eqn{a \geq 1}, \eqn{b \geq 1}.
}

\value{
  The family returns an object of class \code{"gamlss2.family"}.

  Function \code{shiftlog()} returns a link specification
  object of class \code{"link-glm"}.
}

\references{
  Kumaraswamy P (1980).
  A Generalized Probability Density Function for Double-Bounded Random Processes.
  \emph{Journal of Hydrology}, 46(1), 79-88.
  \doi{https://doi.org/10.1016/0022-1694(80)90036-0}
}

\seealso{
  \code{\link{gamlss2.family}}, \code{\link{gamlss2}}
}

\examples{## create family object with
## different link specifications
fam <- Kumaraswamy(a.link = shiftlog, b.link = "log")

## simulate data
set.seed(123)
n <- 1000
d <- data.frame("x" = runif(n, -pi, pi))

## true parameters
par <- data.frame(
  "a" = exp(1.2 + sin(d$x)) + 1,
  "b" = 1
)

## sample response
d$y <- fam$r(1, par)

## estimate model using the Kumaraswamy family
## with shiftlog and log-link
b <- gamlss2(y ~ s(x), data = d,
  family = Kumaraswamy(a.link = shiftlog, b.link = "log"))

## plot estimated effect
plot(b)

## plot residual diagnostics
plot(b, which = "resid")
}

\keyword{regression}
\keyword{models}
\keyword{distribution}
