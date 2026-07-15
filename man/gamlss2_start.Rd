\name{gamlss2_start}
\alias{gamlss2_start}
\title{Model Initialization and Parameter Constraints}

\description{
  Guidance on supplying starting values through argument \code{start} in
  \code{\link{gamlss2}} and on fixing distribution parameters via
  \code{fixed} in \code{\link{gamlss2_control}}.
}

\details{
  The \code{start} argument can be used to initialize estimation in several ways.
  Supported inputs include:
  \itemize{
    \item a named numeric vector or list giving one starting value per
      distribution parameter, such as \code{mu}, \code{sigma}, \code{nu}, or
      \code{tau};
    \item a data frame or matrix with one column per distribution parameter and
      one row per observation, providing starting values for the full additive
      predictors;
    \item a fitted \code{"gamlss2"} object, in which case \code{fitted()}
      values are used as starting values for all distribution parameters;
    \item a \code{coef.gamlss2} object, allowing coefficient-wise
      initialization of linear and smooth terms.
  }

  Parameters can be fixed during estimation using the \code{fixed} control
  argument. This is useful, for example, when shape parameters should remain
  constant while the remaining parameters are estimated.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2_control}}, \code{\link{RS}}
}

\examples{\dontrun{## distributional model
f <- Ozone ~ s(Wind) | s(Wind) | 1 | 1

## scalar starting values for the distribution parameters
m1 <- gamlss2(f, data = airquality, family = BCT,
  start = c(mu = 200, sigma = 0.1, nu = 0, tau = 10))

## fix shape parameters during estimation
m2 <- gamlss2(f, data = airquality, family = BCT,
  start = c(mu = 200, sigma = 0.1, nu = 0, tau = 10),
  fixed = c(nu = TRUE, tau = TRUE))

## coefficient estimates of the fitted model
coef(m2)

## use fitted additive predictors as starting values
m3 <- gamlss2(f, data = airquality, family = BCT,
  start = fitted(m2))

## same using the fitted model object directly
m4 <- gamlss2(f, data = airquality, family = BCT,
  start = m2)

## coefficient-wise initialization
m5 <- gamlss2(f, data = airquality, family = BCT,
  start = coef(m2))
}}

\keyword{models}
\keyword{regression}
