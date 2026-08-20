\name{Kumaraswamy}
\alias{Kumaraswamy}
\alias{KS}
\alias{shiftlog}

\title{Kumaraswamy Distribution}

\description{
  The two-parameter Kumaraswamy family for continuous responses strictly
  between 0 and 1, together with a shifted logarithmic link specification.
}

\usage{
## the Kumaraswamy family
Kumaraswamy(a.link = shiftlog, b.link = shiftlog, ...)
KS(a.link = shiftlog, b.link = shiftlog, ...)

## the log(mu - shift) link specification
shiftlog(shift = 1)
}

\arguments{
  \item{a.link}{A character string, link object, or function specifying the
    link for shape parameter \code{a}.}
  \item{b.link}{A character string, link object, or function specifying the
    link for shape parameter \code{b}.}
  \item{shift}{A single finite number specifying the lower bound imposed by
    the shifted logarithmic link.}
  \item{\dots}{Currently not used.}
}

\details{
  The Kumaraswamy distribution is continuous on \eqn{(0, 1)} with strictly
  positive shape parameters \eqn{a} and \eqn{b}. Its probability density
  function is
  \deqn{f(y; a, b) = a b y^{a - 1} (1 - y^a)^{b - 1},}
  and its cumulative distribution function is
  \deqn{F(y; a, b) = 1 - (1 - y^a)^b.}
  Consequently, the quantile function is
  \deqn{F^{-1}(p; a, b) =
    \left(1 - (1 - p)^{1/b}\right)^{1/a}.}

  For a parameter \eqn{\mu}, \code{shiftlog(shift = s)} defines the link
  and inverse link
  \deqn{g(\mu) = \log(\mu - s), \qquad
    g^{-1}(\eta) = \exp(\eta) + s,}
  with domain \eqn{\mu > s}. The default \code{shift = 1} therefore
  restricts both shape parameters to be greater than 1. In that case the
  distribution has the unique interior mode
  \deqn{\left(\frac{a - 1}{a b - 1}\right)^{1/a}.}
  When alternative links allow \eqn{a \leq 1} or \eqn{b \leq 1}, the
  family returns \code{NA} for the interior mode.
}

\value{
  \code{Kumaraswamy()} and \code{KS()} return the same
  \code{"gamlss2.family"} object. It provides probability density,
  cumulative distribution, quantile, and random-generation functions, as
  well as the mean, variance, and interior mode. Its random generator follows
  the family convention \code{random(par, n)}.

  \code{shiftlog()} returns a link specification of class \code{"link-glm"}.
}

\references{
  Kumaraswamy P (1980). \dQuote{A Generalized Probability Density Function
  for Double-Bounded Random Processes.}
  \emph{Journal of Hydrology}, 46(1), 79--88.
  \doi{10.1016/0022-1694(80)90036-0}
}

\seealso{
  \code{\link{gamlss2.family}}, \code{\link{gamlss2}}
}

\examples{
## create family object with different link specifications
fam <- Kumaraswamy(a.link = shiftlog, b.link = "log")

## simulate data
set.seed(1328)
n <- 1000
d <- data.frame("x" = sort(runif(n, -pi, pi)))

## true parameters
par <- data.frame(
  "a" = exp(1.2 + sin(d$x)) + 1,
  "b" = 1
)

## sample one response for every parameter row
d$y <- fam$random(par, 1)

## estimate model using the Kumaraswamy family
m <- gamlss2(y ~ s(x), data = d, family = fam)

## plot estimated effect
plot(m)

## plot residual diagnostics
plot(m, which = "resid")

## predict quantiles
p <- quantile(m, probs = c(0.05, 0.1, 0.5))

## plot
plot(d, pch = 19, col = adjustcolor(1, 0.3))
matlines(d$x, p, lty = 1, col = 4, lwd = 2)
axis(4, at = p[1, ], labels = names(p))
}

\keyword{regression}
\keyword{models}
\keyword{distribution}
