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
      predictors on the link scale;
    \item a fitted \code{"gamlss2"} object, in which case \code{fitted()}
      values are used as starting values for all distribution parameters;
    \item a \code{coef.gamlss2} object, allowing coefficient-wise
      initialization of linear and smooth terms.
  }

  Parameters can be fixed during estimation using the \code{fixed} control
  argument. This is useful, for example, when shape parameters should remain
  constant while the remaining parameters are estimated.

  Starting values from a previous fit can reduce optimizer iterations.
  For coefficient-wise initialization, \code{coef(object, full = TRUE)} includes
  smooth coefficients; \code{lambdas = TRUE} also includes smoothing parameters.
  Compare iteration counts using the same model, observations, and controls.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2_control}}, \code{\link{RS}}
}

\examples{\dontrun{air <- subset(airquality, !is.na(Ozone))

## hold the smoothing parameter fixed for comparable refits
f <- Ozone ~ s(Wind, sp = 1) | 1
ctrl <- gamlss2_control(maxit = 100, eps = 1e-7)

## scalar starting values for the distribution parameters
m1 <- gamlss2(f, data = air, family = GA,
  start = c(mu = 200, sigma = 0.1), control = ctrl)

## reuse fitted additive predictors on the link scale
m2 <- gamlss2(f, data = air, family = GA,
  start = fitted(m1, type = "link"), control = ctrl)

## coef() without full = TRUE initializes only the linear coefficients
m3 <- gamlss2(f, data = air, family = GA,
  start = coef(m1), control = ctrl)

## reuse all linear and smooth coefficients from the previous fit
m4 <- gamlss2(f, data = air, family = GA,
  start = coef(m1, full = TRUE, lambdas = TRUE, dropall = FALSE),
  control = ctrl)

## full coefficient starts reduce iterations in this example
## check convergence and compare the final log-likelihoods as well
fits <- list(scalar = m1, predictors = m2, linear = m3, full = m4)
data.frame(
  iterations = sapply(fits, function(m) m$iterations),
  converged = sapply(fits, function(m) m$converged),
  logLik = sapply(fits, function(m) as.numeric(logLik(m)))
)

## fix BCT shape parameters at their supplied starting values
m5 <- gamlss2(Ozone ~ s(Wind, sp = 1) | 1 | 1 | 1,
  data = air, family = BCT,
  start = c(mu = 200, sigma = 0.1, nu = 1, tau = 10),
  fixed = c(nu = TRUE, tau = TRUE), maxit = 100)

## fixed flags must also be supplied when reusing coefficients
m6 <- gamlss2(Ozone ~ s(Wind, sp = 1) | 1 | 1 | 1,
  data = air, family = BCT,
  start = coef(m5, full = TRUE, dropall = FALSE),
  fixed = c(nu = TRUE, tau = TRUE), maxit = 100)
unique(fitted(m6, parameter = c("nu", "tau")))
}}

\keyword{models}
\keyword{regression}

