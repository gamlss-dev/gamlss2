\name{gamlss2} 
\alias{gamlss2}
\alias{gamlss2.formula}
\alias{gamlss2.list}

\title{Generalized Additive Models for Location Scale and Shape}

\description{
  Estimation of generalized additive models for location scale and shape (GAMLSS).
  The model fitting function \code{gamlss2()} provides flexible infrastructures
  to estimate the parameters of a response distribution. The number of distributional
  parameters is not fixed, see \code{\link{gamlss2.family}}. Moreover, \code{gamlss2()}
  supports all smooth term constructors from the \pkg{mgcv} package in addition to the
  classical model terms as provided by \pkg{gamlss} and \code{gamlss.add}.
}

\usage{
gamlss2(x, \dots)

\method{gamlss2}{formula}(formula, data, family = NO,
  subset, na.action, weights, offset, start = NULL,
  control = gamlss2_control(\dots), \dots)

\method{gamlss2}{list}(x, \dots)
}

\arguments{
  \item{formula}{A GAM-type \code{\link[stats]{formula}} or \code{\link[Formula]{Formula}}. All
    smooth terms of the \pkg{mgcv} package are supported, see also \code{\link[mgcv]{formula.gam}}.}
  \item{x}{For \code{gamlss.list()} \code{x} is a list of \code{\link[stats]{formula}}s.}
  \item{data}{A data frame or list or environment containing the variables in the model.
    If not found in \code{data}, the variables are taken from \code{environment(formula)},
    typically the environment from which \code{gamlss2} is called.}
  \item{family}{A \code{gamlss.family} or \code{gamlss2.family} object used to define
    distribution and the link functions of the parameters.}
  \item{subset}{An optional vector specifying a subset of observations to be used
    in the fitting process.}
  \item{na.action}{\code{NA} processing for setting up the \code{\link[stats]{model.frame}}.}
  \item{weights}{An optional vector of prior \code{weights} to be used in the fitting process.
    Should be \code{NULL} or a numeric vector.}
  \item{offset}{This can be used to specify an \emph{a priori} known components to be
    included in the linear predictors during fitting. Please note that if only a single
    numeric vector is provided, the offset will be assigned to the first specified parameter
    of the distribution. In the case of multiple offsets, a data frame or list must be supplied.
    Each offset is assigned in the same order as the parameters of the distribution
    specified in the family object.}
  \item{start}{Starting values for estimation algorithms.}
  \item{control}{A list of control arguments, see \code{\link{gamlss2_control}}.}
  \item{\dots}{Arguments passed to \code{\link{gamlss2_control}}.}
}

\details{
  The model fitting function \code{gamlss2()} provides flexible infrastructures for the
  estimation of GAMLSS.
  \itemize{
    \item Distributional models are specified using family objects, either from the
      \pkg{gamlss.dist} package or using \code{\link{gamlss2.family}} objects.
    \item Estimation is carried out through a Newton-Raphson/Fisher scoring algorithm,
      see function \code{\link{RS}}. The estimation algorithms can also be exchanged using 
      \code{\link{gamlss2_control}}. Additionally, if an \code{optimizer} is specified by the
      family object, this optimizer function will be employed for estimation.
    \item The return value is determined by the object returned from the optimizer function,
      typically an object of class \code{"gamlss2"}. Default methods and extractor functions
      are available for this class. Nevertheless, users have the flexibility to supply
      their own optimizer function, along with user-specific methods tailored for the
      returned object.
  }
} 

\value{
  The return value is determined by the object returned from the optimizer function.
  By default, the optimization is performed using the \code{\link{RS}} optimizer function
  (see \code{\link{gamlss2_control}}), yielding an object of class \code{"gamlss2"}.
  Default methods and extractor functions are available for this class.
}

\references{
  Rigby, R.A. and Stasinopoulos, D.M. (2005). Generalized Additive Models for Location,
  Scale and Shape (with Discussion). \emph{Journal of the Royal Statistical Society},
  Series C (Applied Statistics), \bold{54}: 507--554.
  \doi{10.1111/j.1467-9876.2005.00510.x}

  Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019).
	\emph{Distributions for Modeling Location, Scale, and Shape: Using GAMLSS in R},
  Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.
  \doi{10.1201/9780429298547}

  Stasinopoulos, D. M., and Rigby, R. A. (2007). Generalized Additive Models for
  Location Scale and Shape (GAMLSS) in R. \emph{Journal of Statistical Software},
  \bold{23}(7), 1--46. \doi{doi.org/10.18637/jss.v023.i07}

  Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani, F. (2017).
  \emph{Flexible Regression and Smoothing: Using GAMLSS in R},
  Chapman and Hall/CRC. \doi{10.1201/b21973}

  (see also \url{https://www.gamlss.com/}).
}

\seealso{
  \code{\link{RS}}, \code{\link{gamlss2_control}}, \code{\link{gamlss2.family}}
}

\examples{## Load the abdominal circumference data
data("abdom", package = "gamlss.data")

## Specify the model Formula.
f <- y ~ s(x) | s(x) | s(x) | s(x)

## Estimate model.
b <- gamlss2(f, data = abdom, family = BCT)

## Model summary.
summary(b)

## Plot estimated effects.
plot(b, which = "effects")

## Plot diagnostics.
plot(b, which = "resid")

## Predict parameters.
par <- predict(b)

## Predict quantiles.
pq <- NULL
for(q in c(0.05, 0.5, 0.95))
  pq <- cbind(pq, family(b)$q(q, par))

## Plot.
plot(y ~ x, data = abdom, pch = 19,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
matplot(abdom$x, pq, type = "l", lwd = 2,
  lty = 1, col = 4, add = TRUE)

## Use of starting values.
m <- gamlss2(f, data = abdom, family = BCT,
  start = c(mu = 200, sigma = 0.1, nu = 0, tau = 10))

## Fix some parameters.
m <- gamlss2(f, data = abdom, family = BCT,
  start = c(mu = 200, sigma = 0.1, nu = 0, tau = 10),
  fixed = c(nu = TRUE, tau = TRUE))

## Estimated coefficients (intercepts).
coef(m)

## Starting values using full predictors.
m <- gamlss2(f, data = abdom, family = BCT,
  start = fitted(m))

## Same with.
m <- gamlss2(f, data = abdom, family = BCT,
  start = m)
}

\keyword{models}
\keyword{regression}

