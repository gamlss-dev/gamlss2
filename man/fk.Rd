\name{fk}
\alias{fk}
\alias{special_fit.fk}
\alias{special_predict.fk.fitted}

\title{Free-Knot Regression Splines for GAMLSS}

\description{
  Constructs a univariate regression spline whose interior knot locations,
  and optionally the number of interior knots, are selected from the current
  weighted working response.
}

\usage{
fk(x, k = 10L, nk = NULL, degree = 1L,
  criterion = "bic", knots = NULL, K = NULL,
  reselect = FALSE, optim.control = list(), ...)
}

\arguments{
  \item{x}{
    A numeric covariate.
  }
  \item{k}{
    Integer, the largest number of free interior knots
    considered when \code{nk = NULL}. Candidate models from zero through
    \code{k} interior knots are compared.
  }
  \item{nk}{
    Integer, optional fixing the number of free interior knots.
    The locations are still estimated.
  }
  \item{degree}{
    Integer, the polynomial degree of the B-spline basis.
    The default \code{1} gives a piecewise-linear spline.
  }
  \item{criterion}{
    Character selection criterion. Available choices are \code{"bic"},
    \code{"aic"}, \code{"gaic"}, \code{"aicc"}, and \code{"gcv"}.
  }
  \item{knots}{
    Optional numeric vector of fixed interior knot locations. The values must
    be distinct and strictly inside the range of \code{x}. When \code{x}
    is transformed, locations are specified on the transformed scale.
    Supplying \code{knots} overrides automatic location optimization.
  }
  \item{K}{
    Optional nonnegative penalty multiplier used by
    \code{criterion = "gaic"}. The model-level value, or \code{2}, is used
    when \code{NULL}.
  }
  \item{reselect}{
    Logical. If \code{FALSE}, the number of knots selected during the first
    backfitting update is retained and its locations are warm-started in
    later updates. If \code{TRUE}, the number is selected again in every
    update.
  }
  \item{optim.control}{
    Named list passed as the \code{control} argument to
    \code{\link[stats]{optim}} for knot-location optimization.
  }
  \item{\dots}{
    Additional entries merged into \code{optim.control}.
  }
}

\details{
  For each candidate number of interior knots, \code{fk()} uses a B-spline
  basis with fixed, repeated boundary knots. Interior knots are represented
  by positive ordered gaps. This removes permutation ambiguity and prevents
  knot crossings during optimization.

  Conditional on the knot locations, spline coefficients are profiled out
  by weighted least squares. Knot optimization therefore only operates on
  the interior locations. The fitted model term is centered.

  When \code{nk = NULL}, models with zero through \code{k} interior knots
  are compared using \code{criterion}. Effective degrees of freedom include
  both the centered spline coefficients and the estimated interior knot
  locations. Fixed locations supplied through \code{knots} do not add
  knot-location degrees of freedom.
}

\value{
  The constructor returns an object of classes \code{"special"} and
  \code{"fk"}.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{special_terms}}
}

\examples{\dontrun{## load data
data("SpirometryUS")

## subset for female
d <- subset(SpirometryUS, gender == "female")

## Automatically select zero through four free interior knots
m1 <- gamlss2(fev1 ~ fk(age, k = 4) | ., data = d, family = NO)

## estimated effects
plot(m1)

## predict quantiles
qu <- c(0.025, seq(0.1, 0.9, by = 0.1), 0.975)
fit <- quantile(m1, probs = qu)

## plot
plot(fev1 ~ age, data = d)
i <- order(d$age)
matlines(d$age[i], fit[i, ],
  lty = 1, lwd = 2, col = c(2, rep(4, ncol(fit) - 1)))
}}

\keyword{models}
\keyword{smooth}

