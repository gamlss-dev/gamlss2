\name{gamlss2.family} 
\alias{gamlss2.family}
\alias{family.gamlss2}

\title{Family Objects in \pkg{gamlss2}}

\description{
  A \code{"gamlss2.family"} object describes a response distribution for use with
  \code{\link{gamlss2}}. It tells the fitting algorithm which distributional
  parameters exist, which link functions they use, and which distributional
  functions are available, such as densities, quantiles, random-number generators,
  or score functions.
}

\usage{
\method{family}{gamlss2}(object, ...)
}

\arguments{
  \item{object}{A fitted \code{"gamlss2"} model object.}
  \item{\dots}{Not used currently.}
}

\details{
  A \pkg{gamlss2} family object is a list of class \code{"gamlss2.family"}.
  It defines the response distribution and the functions needed to fit and use
  that distribution in a model.

  The minimum required elements are:

  \itemize{
    \item \code{family}: the name of the distribution (character string).
    \item \code{names}: a character vector of parameter names (e.g., \code{c("mu", "sigma")}).
    \item \code{links}: a named character vector specifying the link function for each parameter
          (e.g., \code{c(mu = "identity", sigma = "log")}), or a list of link functions, see
          \code{\link{softplus}}.
    \item \code{pdf(y, par, log = FALSE, ...)}: a function to evaluate the (log-)density.
  }

  The \code{pdf()} function must accept the response \code{y}, a named list
  \code{par} of evaluated parameter values (e.g., \code{par$mu}, \code{par$sigma}),
  a logical \code{log}, and optional additional arguments.

  Further elements can be added depending on how fully the family should support
  fitting, prediction, diagnostics, and simulation. Common optional elements are:

  \itemize{
    \item \code{score}: a named list of functions (one per parameter), each computing the
      first derivative of the log-likelihood with respect to the linear predictor:
      \code{score[[param]](y, par, ...)}.
    \item \code{hess}: a named list of functions computing second derivatives (the negative Hessian).
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
    \item \code{initialize}: a named list of initialization functions, one for each parameter (e.g.,
      \code{initialize$mu(y, ...)}), used to generate starting values.
    \item \code{valid.response(x)}: a function that checks whether the response is valid (e.g.,
      numeric, non-factor).
    \item \code{optimizer()}: an optional function defining a custom optimization
      method for use with \code{\link{gamlss2}}, see also \code{\link{RS}}.
  }

  If analytical score or Hessian functions are not provided, they are approximated
  numerically. If quantile residuals are required, a \code{cdf()} function must
  be available.
}

\value{
  For \code{family.gamlss2()}, the family object stored in a fitted
  \code{"gamlss2"} model is returned. A \code{"gamlss2.family"} object itself is
  a list describing the response distribution, its parameters, their link
  functions, and the distribution-specific methods used during fitting,
  prediction, diagnostics, and simulation.
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

## specify the model formula
f <- y ~ s(x) | s(x)

## estimate model
m1 <- gamlss2(f, data = abdom, family = Normal)

## plot estimated effects
plot(m1, which = "effects")

## plot diagnostics
plot(m1, which = "resid")

## predict parameters
par <- predict(m1)

## predict quantiles
pq <- sapply(c(0.05, 0.5, 0.95), function(q) family(m1)$quantile(q, par))

## visualize
plot(y ~ x, data = abdom, pch = 19,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
matplot(abdom$x, pq, type = "l", lwd = 2,
  lty = 1, col = 4, add = TRUE)

## another example that defines only the density
## function; in this case all derivatives are
## approximated numerically. For residual diagnostics,
## the $cdf() and $quantile() functions are needed as well
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
m2 <- gamlss2(f, data = rent, family = Gamma)

## visualize estimated effects
plot(m2, which = "effects")

## diagnostics, needs the $cdf() and $quantile() function!
plot(m2, which = "resid")
}}

\keyword{regression}
\keyword{models}
\keyword{distribution}
