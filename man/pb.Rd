\name{pb}
\alias{pb}

\title{P-Splines for GAMLSS}

\description{
  Estimation of P-splines using an efficient local maximum likelihood approach to automatically
  select the smoothing parameter. According to the inventors of P-splines,
  \code{pb} stands for "penalized beta" splines or "Paul and Brian".
}

\usage{
pb(x, k = 20, ...)
}

\arguments{
  \item{x}{The variable that should be used for estimation.}
  \item{k}{The dimension of the B-spline basis to represent the smooth term.}
  \item{\dots}{Further arguments passed to function \code{\link[mgcv]{s}}.}
}

\value{
  The function returns a smooth specification object of class \code{"ps.smooth.spec"}, see
  also \code{\link[mgcv]{smooth.construct.ps.smooth.spec}}.
}

\details{
  Function \code{pb()} is an \bold{internal} wrapper function that calls \code{\link[mgcv]{s}}
  to set up a smooth specification object that can be used for model fitting with
  \code{\link{gamlss2}}. Using \code{pb()}, an efficient local maximum likelihood
  approach is used to estimate the smoothing parameter. See the reference for details.
}

\references{
  Eilers PHC, Marx BD (1996). \dQuote{Flexible Smoothing with B-Splines and Penalties.}
  \emph{Statistical Science}, \bold{11}(2), 89--121.
  \doi{10.1214/ss/1038425655}

  Rigby RA, Stasinopoulos DM (2014).
  \dQuote{Automatic Smoothing Parameter Selection in GAMLSS with an Application to Centile Estimation.}
  \emph{Statistical Methods in Medical Research}, \bold{23}(4),  318--332.
  \doi{10.1177/0962280212473302}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link[mgcv]{smooth.construct.ps.smooth.spec}}
}

\examples{
\dontshow{ if(!requireNamespace("gamlss.data")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }}

## Load head circumference data.
data("dbhh", package = "gamlss.data")

## Specify the model Formula.
f <- head ~ pb(age) | pb(age) | pb(age) | pb(age)

## Estimate model.
b <- gamlss2(f, data = dbhh, family = BCT)

## Plot estimated effects.
plot(b, which = "effects")

## Plot diagnostics.
plot(b, which = "resid")

## Predict parameters.
par <- predict(b)

## Predict quantiles.
pq <- NULL
for(q in seq(0.05, 0.95, by = 0.1))
  pq <- cbind(pq, family(b)$q(q, par))

## Plot.
plot(head ~ age, data = dbhh, pch = 19,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
matplot(dbhh$age, pq, type = "l",
  lty = 1, col = 4, add = TRUE)
}

\keyword{utilities}
\keyword{regression}

