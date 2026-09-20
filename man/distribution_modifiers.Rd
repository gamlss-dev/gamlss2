\name{distribution_modifiers}
\alias{distribution_modifiers}
\alias{transform_family}
\alias{logit_family}
\alias{truncate_family}
\alias{censor_family}
\alias{inflate_family}
\alias{hurdle_family}
\alias{discretize_family}
\alias{round_family}

\title{Construct Modified Distribution Families}

\description{
  Constructs \code{"gamlss2.family"} objects by transforming, conditioning,
  coarsening, or adding point masses to an existing distribution family.
  Modified families can be used directly in \code{\link{gamlss2}} models and
  can be nested to compose several operations.
}

\usage{
transform_family(family = NO, transform, inverse, log_jacobian,
  increasing = TRUE, support = c(-Inf, Inf), valid.response = NULL,
  name = NULL)

logit_family(family = NO)

truncate_family(family = NO, lower = -Inf, upper = Inf,
  include = c(FALSE, TRUE))

censor_family(family = NO, lower = -Inf, upper = Inf)

inflate_family(family = NO, at = 0)

hurdle_family(family = PO, at = 0)

discretize_family(family = NO, breaks,
  values = seq_len(length(breaks) - 1L), right = TRUE)

round_family(family = NO, digits = 0L)
}

\arguments{
  \item{family}{
    A distribution family accepted by the \pkg{gamlss2} family
    infrastructure, such as \code{NO} or \code{PO}, or another modified
    family.
  }
  \item{transform}{
    A fixed monotone function mapping observations from the base scale to
    the response scale.
  }
  \item{inverse}{
    The inverse function of \code{transform}.
  }
  \item{log_jacobian}{
    A function returning the logarithm of the absolute derivative of
    \code{inverse}, evaluated on the response scale.
  }
  \item{increasing}{
    Logical. Is \code{transform} increasing?
  }
  \item{support}{
    Numeric vector with the two response-scale support bounds. These bounds
    are also used when evaluating the CDF outside the transformed support.
  }
  \item{valid.response}{
    Optional function returning logical values that indicate whether
    observations lie in the transformed support.
  }
  \item{name}{
    Optional character string used as the modified family label.
  }
  \item{lower, upper}{
    Fixed scalar lower and upper truncation or censoring boundaries.
  }
  \item{include}{
    Logical vector of length two indicating whether the lower and upper
    truncation boundaries are retained. The default interval is
    \eqn{(lower, upper]}.
  }
  \item{at}{
    One or more fixed inflation points, or one fixed hurdle point.
  }
  \item{breaks}{
    Strictly increasing numeric bin boundaries beginning at \code{-Inf} and
    ending at \code{Inf}.
  }
  \item{values}{
    Strictly increasing, finite values assigned to the bins.
  }
  \item{right}{
    Logical. Should bins be right-closed intervals \eqn{(a, b]} rather than
    left-closed intervals \eqn{[a, b)}?
  }
  \item{digits}{
    Integer number of decimal places retained by rounding.
  }
}

\details{
  For a fixed monotone transformation \eqn{Y = g(X)},
  \code{transform_family} uses
  \deqn{\log f_Y(y) = \log f_X(g^{-1}(y)) + \log |d g^{-1}(y) / dy|.}
  Base-family parameter scores, Hessians, and optimized updates can therefore
  be evaluated directly at the inverse-transformed response. Increasing and
  decreasing transformations are supported. The convenience function
  \code{logit_family} constructs the logistic push-forward,
  which has strict support \eqn{0 < Y < 1}.
  Exact zeros and ones are not valid observations for this
  continuous transformation.

  The function \code{truncate_family} conditions the base distribution on
  the interval selected by \code{lower}, \code{upper}, and \code{include}.
  Its density or probability mass is divided by the probability of the
  retained interval. In contrast, \code{censor_family} observes
  \deqn{Y = \min(\max(X, lower), upper)}
  and consequently creates point masses at finite censoring boundaries.
  Truncation and censoring reuse base-family derivatives where possible and
  differentiate only the additional probability terms.

  For inflation points \eqn{c_j}, \code{inflate_family} constructs
  \deqn{F(y) = \pi_0 F_0(y) + \sum_j \pi_j I(c_j \le y).}
  Natural base-distribution mass at an inflation point is retained. The
  appended parameters \code{pi1}, \code{pi2}, and so on have log links and
  represent positive odds \eqn{w_j} relative to the base component. Thus,
  \deqn{\pi_0 = \frac{1}{1 + \sum_j w_j}, \qquad \pi_j = \frac{w_j}{1 + \sum_\ell w_\ell}.}
  This is the baseline-odds parameterization also used by
  \code{\link{mixture_family}}. Quantiles account for every probability jump.

  The function \code{hurdle_family} instead constructs
  \deqn{F(y) = \pi_1 I(c \le y) + \pi_0 F_0(y \mid X \ne c).}
  Hence, the hurdle probability is the complete mass at \eqn{c}, including
  when the base distribution already assigns positive mass to that point.

  The function \code{discretize_family} assigns the probability of every
  latent interval defined by \code{breaks} to the corresponding element of
  \code{values}. The function \code{round_family} uses intervals centered on
  a decimal grid and is restricted to continuous base distributions. Its
  quantile and random-number functions follow R's \code{round()} convention.

  Inflated, hurdle, and fixed-bin families provide exact means and variances
  when the required base-family moments are available. Generic censored and
  rounded moments require partial moments, which are not part of the family
  interface, and therefore report an informative error rather than applying
  an incorrect numerical approximation.

  Modified continuous distributions containing explicit point masses are
  treated as mixed distributions and provide randomized quantile residuals.
  Bounds, atom locations, breaks, and transformations are fixed when the
  family is constructed; they are not modeled distributional parameters.

  Arbitrary weighted distributions, splices, non-monotone transformations,
  order-statistic distributions, convolutions, and compound distributions
  are not provided. Generic implementations would require numerical
  normalization or derivatives inside fitting iterations and could not be
  supplied reliably by the existing family interface.
}

\value{
  An object of class \code{"gamlss2.family"}. Existing base-family parameters
  retain their names, order, and links. Inflation and hurdle families append
  positive-odds parameters named \code{pi1}, \code{pi2}, and so on.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2.family}},
  \code{\link{mixture_family}}, \code{\link{distributions3_family}}
}

\examples{\dontrun{## simulate data
set.seed(6020)
n <- 300
x <- sort(runif(n, -3, 3))

## logistic transformation of a normal distribution
fam <- logit_family(NO)
y <- fam$random(list(mu = sin(x), sigma = exp(-1 + cos(x))), n)
hist(y)
m <- gamlss2(y ~ s(x) | s(x), family = fam)
plot(m)

## zero-truncated Poisson distribution
fam <- truncate_family(PO, lower = 0)
y <- fam$random(list(mu = exp(1/2 + sin(x))), n)
m <- gamlss2(y ~ s(x), family = fam)
plot(x, y)

## small helper for plotting fitted quantiles
fq <- function(model) {
  matlines(x,
    quantile(model, probs = c(0.05, 0.5, 0.95)),
    lty = 1, lwd = 2, col = 4
  )
}

fq(m)

## normal distribution censored at both boundaries
fam <- censor_family(NO, lower = -1, upper = 1)
y <- fam$random(list(mu = sin(x), sigma = exp(-1 + cos(x))), n)
hist(y)
m <- gamlss2(y ~ s(x) | s(x), family = fam)
plot(m)
plot(x, y)
fq(m)

## zero-inflated Poisson distribution
fam <- inflate_family(PO, at = 0)
y <- fam$random(list(mu = exp(1/2 + sin(x)), pi1 = plogis(cos(x))), n)
barplot(table(y))
m <- gamlss2(y ~ s(x) | s(x), family = fam)
plot(m)
plot(x, y)
fq(m)

## zero-hurdle Poisson distribution
fam <- hurdle_family(PO, at = 0)
y <- fam$random(list(mu = exp(1/2 + sin(x)), pi1 = plogis(cos(x))), n)
barplot(table(y))
m <- gamlss2(y ~ s(x) | s(x), family = fam)
plot(m)
plot(x, y)
fq(m)

## fixed-bin discretization and decimal rounding
fam <- discretize_family(NO, breaks = c(-Inf, -2, -1, 0, 1, 2, Inf))
y <- fam$random(list(mu = sin(x), sigma = exp(-1 + cos(x))), n)
barplot(table(y))
m <- gamlss2(y ~ s(x) | s(x), family = fam)
plot(m)
plot(x, y)
fq(m)

fam <- round_family(NO, digits = 1)
y <- fam$random(list(mu = sin(x), sigma = exp(-1 + cos(x))), n)
m <- gamlss2(y ~ s(x) | s(x), family = fam)
plot(m)
plot(m, which = "resid")
}}

\keyword{distribution}
\keyword{models}
\keyword{regression}

