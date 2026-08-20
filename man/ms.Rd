\name{ms}
\alias{ms}
\alias{special_fit.ms}
\alias{special_predict.ms.fitted}

\title{Monotonic P-Splines for GAMLSS}

\description{
  Constructs a monotonic penalized B-spline (P-spline) model term for use
  with \code{gamlss2}. The smooth can be constrained to be monotonically
  increasing or decreasing. The smoothing parameter is either supplied by
  the user or selected automatically during estimation.
}
\usage{
ms(x, mono = c("up", "down"), k = 10,
  criterion = "aicc", lambda = NULL,
  lambda.min = 1e-8, lambda.max = 1e8,
  lambda.tol = 1e-4, ...)
}
\arguments{
  \item{x}{
    A numeric covariate to be represented by a monotonic P-spline.
  }
  \item{mono}{
    Character specifying the direction of monotonicity. Use
    \code{"up"} for a monotonically increasing effect or \code{"down"}
    for a monotonically decreasing effect.
  }
  \item{k}{
    Integer giving the dimension of the B-spline basis. The default is
    \code{10}.
  }
  \item{criterion}{
    Character specifying the criterion used for automatic smoothing
    parameter selection. Available criteria are \code{"aicc"},
    \code{"aic"}, \code{"gaic"}, \code{"bic"}, and \code{"gcv"}.
    The default is \code{"aicc"}.
  }
  \item{lambda}{
    Optional positive smoothing parameter. If \code{NULL}, the smoothing
    parameter is selected automatically. If supplied, \code{lambda} is
    kept fixed during fitting.
  }
  \item{lambda.min, lambda.max}{
    Positive lower and upper limits for automatic smoothing parameter
    selection. The search is performed on the logarithmic scale.
  }
  \item{lambda.tol}{
    Tolerance used in the one-dimensional optimization of the smoothing
    parameter.
  }
  \item{\dots}{
    Further arguments passed to \code{\link[mgcv]{s}} when constructing
    the P-spline basis.
  }
}

\details{
  Function \code{ms()} defines a special model term for \code{gamlss2}.
  The underlying basis is an \code{mgcv} P-spline basis constructed with
  \code{bs = "ps"}.

  For an increasing smooth, adjacent B-spline coefficients are constrained
  to satisfy
  \deqn{\beta_{r+1} - \beta_r \ge 0,}
  whereas for a decreasing smooth the inequalities are reversed. Because
  B-spline basis functions are non-negative and form a partition of unity,
  ordering the spline coefficients in this way produces a monotonic fitted
  effect.

  Estimation is based on penalized weighted least squares. For a given
  smoothing parameter \eqn{\lambda}, the coefficient vector minimizes
  \deqn{
    (z - X\beta)^\top W (z - X\beta)
    + \lambda \beta^\top S \beta
  }
  subject to the monotonicity constraints. The constrained problem is solved
  with an active-set algorithm. Weighted cross-products are computed once
  per backfitting update, so smoothing parameter selection only requires
  operations on matrices having dimension determined by \code{k}.

  If \code{lambda = NULL}, the smoothing parameter is selected by a
  one-dimensional search on \code{log(lambda)} using the criterion specified
  by \code{criterion}.

  The fitted term is centered before it is returned to the additive
  backfitting algorithm.
}

\value{
  The constructor \code{ms()} returns an object of classes
  \code{"special"} and \code{"ms"} containing the P-spline design and
  penalty matrices and information required for fitting and prediction.

  After fitting with \code{gamlss2}, the corresponding fitted special term
  contains the following main components:
  \itemize{
    \item \code{coefficients}: estimated B-spline coefficients.
    \item \code{fitted.values}: centered fitted contribution of the term.
    \item \code{edf}: effective degrees of freedom.
    \item \code{lambdas}: selected or fixed smoothing parameter.
    \item \code{vcov}: covariance/influence matrix conditional on the
          active constraints.
    \item \code{df}: residual degrees of freedom for the term,
          \code{n - edf}.
  }
  Additional components are stored for prediction and for warm starts in
  subsequent backfitting iterations.
}

\references{
  Eilers, P. H. C. and Marx, B. D. (1996).
  \dQuote{Flexible Smoothing with B-Splines and Penalties.}
  \emph{Statistical Science}, 11(2), 89--121.
  \doi{10.1214/ss/1038425655}

  Bollaerts, K., Eilers, P. H. C. and van Mechelen, I. (2006).
  \dQuote{Simple and Multiple P-splines Regression with Shape Constraints.}
  \emph{British Journal of Mathematical and Statistical Psychology},
  59(2), 451--469. \doi{10.1348/00071105X84293}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{special_terms}}
}

\examples{\dontrun{## air quality data
air <- subset(airquality, !is.na(Ozone))
air <- air[order(air$Temp), ]

## monotone Temp effect for mu and sigma
## Wind is modeled by an unconstrained smooth
f <- Ozone ~ s(Wind) + ms(Temp, "up") | ms(Temp, "down")

## estimate model
m <- gamlss2(f, data = air, family = GA)

## model summary
summary(m)

## visualize estimated effects
plot(m)
}}

\keyword{models}
\keyword{smooth}

