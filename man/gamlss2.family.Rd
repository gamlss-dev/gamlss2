\name{gamlss2.family} 
\alias{gamlss2.family}
\alias{family.gamlss2}

\title{Family Objects in \pkg{gamlss2}}

\description{
  Family objects within the package \pkg{gamlss2} are used to specify the
  information required to use a model fitting function. This includes details
  such as parameter names, corresponding link functions, the density function,
  log-likelihood function and derivatives of the log-likelihood with respect to
  the predictors. In addition, family objects are used in the calculation of
  post-modeling statistics, such as residual diagnostics and random number generation.
  An overview can be found in the accompanying details and examples.
}

\details{
  A \pkg{gamlss2} family object must be a \code{list} of class \code{"gamlss2.family"}, containing
  the necessary elements to define the response distribution and how it is handled during model
  fitting. The minimum requirements are as follows:

  \itemize{
    \item \code{family}: the name of the distribution (character string).
    \item \code{names}: a character vector of parameter names (e.g., \code{c("mu", "sigma")}).
    \item \code{links}: a named character vector specifying the link function for each parameter
          (e.g., \code{c(mu = "identity", sigma = "log")}), or a list of link functions, e.g., see
          \code{\link{softplus}}.
    \item \code{pdf(y, par, log = FALSE, ...)}: a function to evaluate the (log-)density.
  }

  The \code{pdf()} function must accept the response \code{y}, a named list \code{par} of evaluated
  parameter values (e.g., \code{par$mu}, \code{par$sigma}), a logical \code{log}, and optional
  additional arguments.

  Optionally, the family object may include:

  \itemize{
    \item \code{score}: a named list of functions (one per parameter), each computing the
      first derivative of the log-likelihood with respect to the linear predictor:
      \code{score[[param]](y, par, ...)}.
    \item \code{hess}: a named list of functions computing second derivatives (negative Hessian).
      For parameters \code{mu} and \code{sigma}, this includes:
      \code{hess[["mu"]]()}, \code{hess[["sigma"]]()}, and optionally
      cross derivatives like \code{hess[["mu.sigma"]]()}.
    \item \code{loglik(y, par, ...)}: a function computing the total log-likelihood.
    \item \code{cdf(y, par, ...)}: cumulative distribution function.
    \item \code{quantile(p, par, ...)}: quantile function.
    \item \code{random(n, par, ...)}: random number generator.
    \item \code{mean(par, ...)}: mean function.
    \item \code{variance(par, ...)}: variance function.
    \item \code{skewness(par, ...)}, \code{kurtosis(par, ...)}: optional higher-order moment functions.
    \item \code{initialize}: a named list of initialization functions, one per parameter (e.g.,
      \code{initialize$mu(y, ...)}), used to generate starting values.
    \item \code{valid.response(x)}: a function that checks whether the response is valid (e.g.,
      numeric, non-factor).
    \item \code{optimizer()}: an optional function to define a custom optimization method for
      use with \code{\link{gamlss2}}, see also \code{link{RS}}.
  }

  If the analytical score or Hessian functions are not provided, they will be approximated
  numerically. If quantile residuals are to be computed, a \code{cdf()} function
  must be provided.
}

\seealso{
  \code{\link{gamlss2}}
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
\dontrun{## new family object for the normal distribution
Normal <- function(...) {
  fam <- list(
    "family" = "Normal",
    "names" = c("mu", "sigma"),
    "links" = c("mu" = "identity", "sigma" = "log"),
    "score" = list(
      "mu" = function(y, par, ...) {
        (y - par$mu) / (par$sigma^2)
      },
      "sigma" = function(y, par, ...) {
        -1 + (y - par$mu)^2 / (par$sigma^2)
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        1 / (par$sigma^2)
      },
      "sigma" = function(y, par, ...) {
        rep(2, length(y))
      },
      "mu.sigma" = function(y, par, ...) {
        rep(0, length(y))
      }
    ),
    "loglik" = function(y, par, ...) {
      sum(dnorm(y, par$mu, par$sigma, log = TRUE))
    },
    "mu" = function(par, ...) {
      par$mu
    },
    "pdf" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "cdf" = function(y, par, ...) {
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    },
    "random" = function(n, par) {
      rnorm(n, mean = par$mu, sd = par$sigma)
    },
    "quantile" = function(p, par) {
      qnorm(p, mean = par$mu, sd = par$sigma)
    },
    "initialize" = list(
      "mu"    = function(y, ...) { (y + mean(y)) / 2 },
      "sigma" = function(y, ...) { rep(sd(y), length(y)) }
    ),
    "mean"      = function(par) par$mu,
    "variance"  = function(par) par$sigma^2,
    "valid.response" = function(x) {
      if(is.factor(x) | is.character(x))
        stop("the response should be numeric!")
      return(TRUE)
    }
  )

  class(fam) <- "gamlss2.family"

  return(fam)
}

## load the abdominal circumference data
data("abdom", package = "gamlss.data")

## specify the model Formula
f <- y ~ s(x) | s(x)

## estimate model
b <- gamlss2(f, data = abdom, family = Normal)

## plot estimated effects
plot(b, which = "effects")

## plot diagnostics
plot(b, which = "resid")

## predict parameters
par <- predict(b)

## predict quantiles
pq <- sapply(c(0.05, 0.5, 0.95), function(q) family(b)$quantile(q, par))

## visualize
plot(y ~ x, data = abdom, pch = 19,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
matplot(abdom$x, pq, type = "l", lwd = 2,
  lty = 1, col = 4, add = TRUE)

## another example using only the density
## function, all derivatives are approximated
## in this case; for residual diagnostics,
## the $cdf() and $quantile() function is needed, too.
Gamma <- function(...) {
  fam <- list(
    "names" = c("mu", "sigma"),
    "links" = c("mu" = "log", "sigma" = "log"),
    "pdf" = function(y, par, log = FALSE, ...) {
      shape <- par$sigma
      scale <- par$mu/par$sigma
      dgamma(y, shape = shape, scale = scale, log = log)
    },
    "cdf" = function(y, par, lower.tail = TRUE, log.p = FALSE) {
      shape <- par$sigma
      scale <- par$mu/par$sigma
      pgamma(y, shape = shape, scale = scale,
        lower.tail = lower.tail, log.p = log.p)
    },
    "quantile" = function(p, par, lower.tail = TRUE, log.p = FALSE) {
      shape <- par$sigma
      scale <- par$mu/par$sigma
       qgamma(p, shape = shape, scale = scale,
         lower.tail = lower.tail, log.p = log.p)
    }
  )

  class(fam) <- "gamlss2.family"

  return(fam)
}

## example using the Munich rent data
data("rent", package = "gamlss.data")

## model formula
f <- R ~ ti(Fl) + ti(A) + ti(Fl, A, bs = "ps") |
  ti(Fl) + ti(A) + ti(Fl, A, bs = "ps")

## estimate model
b <- gamlss2(f, data = rent, family = Gamma)

## visualize estimated effects
plot(b, which = "effects")

## diagnostics, needs the $cdf() and $quantile() function!
plot(b, which = "resid")
}}

\keyword{regression}
\keyword{models}
\keyword{distribution}

