\name{bamlss2}
\alias{bamlss2}

\title{Bayesian GAMLSS Wrapper}

\description{
  A wrapper around \code{\link{gamlss2}} for fitting Bayesian GAMLSS models
  by MCMC sampling. Starting values are computed using \code{maxit}
  backfitting iterations before sampling. The controls \code{n.iter},
  \code{burnin}, and \code{thin} can be specified directly; \code{weights}
  and \code{offset} are passed through to \code{gamlss2}.
}

\usage{
bamlss2(formula, n.iter = 1200, burnin = 200, thin = 1, maxit = 2, \dots)
}

\arguments{
  \item{formula}{A GAM-type \code{\link[stats]{formula}} or
    \code{\link[Formula]{Formula}}. The model can include smooth terms as
    provided by the \pkg{mgcv} package.}
  \item{n.iter}{Integer, the total number of MCMC iterations.}
  \item{burnin}{Integer, the burn-in period.}
  \item{thin}{Integer, the thinning parameter for saved samples.}
  \item{maxit}{Integer, the number of backfitting iterations used to compute
    starting values before running MCMC.}
  \item{\dots}{Further arguments passed to \code{\link{gamlss2}} and/or
    \code{\link{gamlss2_control}}.}
}

\value{
  An object of class \code{"bamlss2"} inheriting from \code{"gamlss2"}
  containing the fitted model and posterior samples.
}

\references{
  Umlauf N, Klein N, Zeileis A (2018).
  \dQuote{BAMLSS: Bayesian Additive Models for Location, Scale, and Shape (and Beyond)}.
  \emph{Journal of Computational and Graphical Statistics}, \bold{27}(3),
  612--627. \doi{10.1080/10618600.2017.1407325}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{mcmc}}, \code{\link{BS}}
}

\examples{\dontrun{
## air quality data
air <- subset(airquality, !is.na(Ozone))

## set the seed
set.seed(123)

## estimate Gamma model with MCMC
m <- bamlss2(
  Ozone ~ s(Temp) + s(Wind) | s(Temp),
  data = air,
  family = GA
)

## posterior summary
summary(m)

## plot estimated effects
plot(m)

## plot samples
plot(m, which = "samples", ask = TRUE)

## predict parameters using samples
pm <- predict(m, type = "parameter", FUN = mean)
print(head(pm))
psd <- predict(m, type = "parameter", FUN = sd)
print(head(psd))
}}

\keyword{models}
\keyword{regression}
