\name{gamlss2.family} 
\alias{gamlss2.family}
\alias{family.gamlss2}

\title{
  Family Objects in \pkg{gamlss2}
}

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
  The following lists the minimum requirements on a \pkg{gamlss2} family object to be used with
  \code{\link{gamlss2}}:
  \itemize{
    \item The family object is expected to return a \code{\link{list}} of class \code{"gamlss2.family"}.
    \item The object must contain the family name as a character string. 
    \item The object must contain the names of the parameters as a character string, as well as
      the corresponding link functions as character string.
    \item The family object must contain a \code{$d()} function to evaluate the (log-)density.
  }
  Furthermore, it is assumed that the density function in a family object has the
  following arguments:

  \code{d(y, par, log = FALSE, ...)}

  where argument \code{y} is the response (possibly a matrix) and \code{par} is a named list
  holding the evaluated parameters of the distribution, e.g., using a normal distribution \code{par}
  has two elements, one for the mean \code{par$mu} and one for the standard deviation
  \code{par$sigma}. The dots argument is for passing special internally used objects, depending
  on the type of model this feature is usually not needed.

  Optionally, the family object holds derivative functions evaluating
  derivatives of the log-likelihood w.r.t. the predictors (or expectations of derivatives).
  For each parameter, these functions must have the following arguments:
 
  \code{function(y, par, ...)}

  for computing the first derivative of the log-likelihood w.r.t. a predictor and

  \code{function(y, par, ...)}

  for computing the _negative_ second derivatives. Within the family object these functions
  are organized in a named list, see the examples below. If these functions are not specified,
  all derivatives will be approximated numerically. Note that also cross derivatives can
  be implemented, e.g., when using the \code{\link{CG}} algorithm for fitting a GAMLSS.

  In addition, for the cumulative distribution function (\code{p(y, par, ...)}), for the quantile
  function (\code{q(y, par, ...)}) or for creating random numbers (\code{r(n, par, ...)}) the same
  structure is assumed.

  Using function \code{\link{gamlss2}} the family objects may
  also specify the \code{optimizer()}er function that
  should be used with this family.
}

\seealso{
  \code{\link{gamlss2}}
}

\examples{\dontrun{## New family object for the normal distribution.
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
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "p" = function(y, par, ...) {
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    },
    "r" = function(n, par) {
      rnorm(n, mean = par$mu, sd = par$sigma)
    },
    "q" = function(p, par) {
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

## Load the abdominal circumference data
data("abdom", package = "gamlss.data")

## Specify the model Formula.
f <- y ~ s(x) | s(x)

## Estimate model.
b <- gamlss2(f, data = abdom, family = Normal)

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

## Another example using only the density
## function, all derivatives are approximated
## in this case. For residual diagnostics,
## the $p() and $q() function is needed, too.
Gamma <- function(...) {
  fam <- list(
    "names" = c("mu", "sigma"),
    "links" = c("mu" = "log", "sigma" = "log"),
    "d" = function(y, par, log = FALSE, ...) {
      shape <- par$sigma
      scale <- par$mu/par$sigma
      dgamma(y, shape = shape, scale = scale, log = log)
    },
    "p" = function(y, par, lower.tail = TRUE, log.p = FALSE) {
      shape <- par$sigma
      scale <- par$mu/par$sigma
      pgamma(y, shape = shape, scale = scale,
        lower.tail = lower.tail, log.p = log.p)
    },
    "q" = function(p, par, lower.tail = TRUE, log.p = FALSE) {
      shape <- par$sigma
      scale <- par$mu/par$sigma
       qgamma(p, shape = shape, scale = scale,
         lower.tail = lower.tail, log.p = log.p)
    }
  )

  class(fam) <- "gamlss2.family"

  return(fam)
}

## Example using the Munich rent data.
data("rent", package = "gamlss.data")

## Model formula.
f <- R ~ ti(Fl) + ti(A) + ti(Fl, A, bs = "ps") |
  ti(Fl) + ti(A) + ti(Fl, A, bs = "ps")

## Estimate model.
b <- gamlss2(f, data = rent, family = Gamma)

## Plot estimated effects.
plot(b, which = "effects")

## Diagnostics, needs the $p() and $q() function!
plot(b, which = "resid")
}}

\keyword{regression}
\keyword{models}
\keyword{distribution}

