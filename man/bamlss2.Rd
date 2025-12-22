\name{bamlss2} 
\alias{bamlss2}

\title{Bayesian GAMLSS Wrapper}

\description{
  The function \code{bamlss2()} is a convenience wrapper around \code{\link{gamlss2}}
  to fit Bayesian GAMLSS models via MCMC sampling. It sets the optimizer to an MCMC
  routine and allows to specify \code{n.iter}, \code{burnin}, and \code{thin} directly.
}

\usage{
bamlss2(formula, n.iter = 1200, burnin = 200, thin = 1, maxit = 2, \dots)
}

\arguments{
  \item{formula}{A GAM-type \code{\link[stats]{formula}} or \code{\link[Formula]{Formula}}.
    The model can include smooth terms as provided by the \pkg{mgcv} package.}
  \item{n.iter}{Integer, the total number of MCMC iterations.}
  \item{burnin}{Integer, the burn-in period.}
  \item{thin}{Integer, thinning parameter for saved samples.}
  \item{maxit}{Integer, number of backfitting iterations to compute starting values
    before running MCMC. If \code{maxit = 0}, sampling starts from default or user-supplied
    starting values.}
  \item{\dots}{Further arguments passed to \code{\link{gamlss2}} and/or
    \code{\link{gamlss2_control}}.}
}

\details{
  Function \code{bamlss2()} calls \code{\link{gamlss2}} and sets the optimizer to a
  Bayesian sampler. By default, a small number of backfitting iterations is performed
  to obtain reasonable starting values for the MCMC sampler.

  The wrapper is designed to behave in line with \code{gamlss()} interfaces, i.e., it
  supports \code{weights} and \code{offset} arguments through \code{\link{gamlss2}}.
}

\value{
  An object of class \code{"bamlss2"} inheriting from \code{"gamlss2"} containing the
  fitted model and posterior samples.
}

\references{
  Umlauf, N., Klein, N., and Zeileis, A. (2018).
  \emph{BAMLSS: Bayesian Additive Models for Location, Scale and Shape (and Beyond)}.
  Journal of Computational and Graphical Statistics, \bold{27}(3), 612--627.
  \doi{10.1080/10618600.2017.1407325}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{mcmc}}, \code{\link{BS}}
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
## load the abdominal circumference data
data("abdom", package = "gamlss.data")

## specify the model Formula
f <- y ~ s(x) | s(x) | s(x) | s(x)

## estimate Bayesian model
b <- bamlss2(f, data = abdom, family = BCT)

## posterior summary
summary(b)

## plot estimated effects
plot(b)

## plot samples
plot(b, which = "samples")

## predict parameters using samples
pm <- predict(b, FUN = mean)
print(head(pm))
psd <- predict(b, FUN = sd)
print(head(psd))
}

