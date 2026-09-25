\name{discretize_family}
\alias{discretize_family}

\title{Discretize a Distribution Family into Fixed Bins}

\description{
  Constructs a discrete distribution family by assigning the probability of
  each fixed interval of a base distribution to a corresponding response
  value.
}

\usage{
discretize_family(family = NO, breaks,
  values = seq_len(length(breaks) - 1L), right = TRUE)
}

\arguments{
  \item{family}{A distribution family accepted by the \pkg{gamlss2} family
    infrastructure.}
  \item{breaks}{A strictly increasing numeric vector of interval boundaries.
    The first and last elements must be \code{-Inf} and \code{Inf}.}
  \item{values}{A strictly increasing vector of finite numeric response values,
    one for each interval defined by \code{breaks}. By default, the intervals
    are labeled consecutively starting at one.}
  \item{right}{Logical. Should the intervals be right-closed, \eqn{(a, b]},
    rather than left-closed, \eqn{[a, b)}?}
}

\details{
  For each interval defined by adjacent elements of \code{breaks},
  \code{discretize_family()} assigns the interval probability under the base
  distribution to the corresponding element of \code{values}. Thus, with
  \code{right = TRUE}, response value \eqn{v_j} has probability
  \deqn{P(X \in (b_j, b_{j+1}])}
  under the base distribution. With \code{right = FALSE}, the intervals are
  left-closed instead.

  The returned family provides probability mass, cumulative distribution,
  quantile, random-number, mean, and variance functions. Its distributional
  parameters and link functions are inherited from the base family.
}

\value{
  An object of class \code{"gamlss2.family"} representing the discretized
  distribution.
}

\seealso{
  \code{\link{family_extensions}}, \code{\link{round_family}},
  \code{\link{gamlss2}}, \code{\link{gamlss2.family}}
}

\examples{\dontrun{
## fixed intervals and their observed values
breaks <- c(-Inf, -1, 0, 1, Inf)
values <- 0:3

## discretize a normal family
fam <- discretize_family(NO, breaks = breaks, values = values)

## simulate data whose location depends on x
set.seed(6020)
n <- 500
x <- sort(runif(n, -3, 3))
par <- data.frame(mu = sin(x), sigma = exp(-0.5 + cos(x)))
y <- fam$random(par, n)
head(y)
d <- data.frame(y = y, x = x)

## plot count data
plot(y ~ x, data = d)

## fitted pdf
fit_family(d$y, family = fam)

## fit a regression model using the discretized family
m <- gamlss2(y ~ s(x) | ., data = d, family = fam)

## inspect the estimated effect and fitted probabilities
plot(m)
d <- predict(m)
head(pdf(d, y))

## predict median
p <- quantile(d, probs = 0.5)
plot(y ~ x, data = d)
lines(p ~ x, lwd = 2, col = 4)
}}

\keyword{distribution}
\keyword{regression}
\keyword{discretization}

