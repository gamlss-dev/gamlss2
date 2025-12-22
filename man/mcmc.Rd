\name{mcmc} 
\alias{mcmc}

\title{Run MCMC Sampling for a Fitted \code{gamlss2} Model}

\description{
  The function \code{mcmc()} runs MCMC sampling for a fitted \code{\link{gamlss2}} model.
  It reuses the fitted model object and merges the MCMC output back into the original
  object such that standard methods (e.g., \code{\link{summary}} and \code{\link{plot}})
  work as expected.
}

\usage{
mcmc(object, n.iter = 1200, burnin = 200, thin = 1)
}

\arguments{
  \item{object}{An object of class \code{"gamlss2"} (or \code{"bamlss2"}) as returned by
    \code{\link{gamlss2}}.}
  \item{n.iter}{Integer, the total number of MCMC iterations.}
  \item{burnin}{Integer, the burn-in period.}
  \item{thin}{Integer, thinning parameter for saved samples.}
}

\details{
  The function sets \code{n.iter}, \code{burnin}, and \code{thin} in the model control
  list and obtains starting values from the fitted object (including smoothing parameters).
  The MCMC sampler \code{\link{BS}} is then called and the resulting samples and posterior
  summaries are merged into the original model object. This approach ensures that the
  result retains \code{terms}, \code{call}, and related components required by extractor
  and plotting functions.

  The function is designed to behave in line with \code{gamlss2()} interfaces with
  respect to \code{weights} and \code{offset} handling. If the fitted model object does
  not contain the required design information for MCMC sampling (e.g., if the model was
  fitted using a reduced \code{control} setting that drops smooth design matrices),
  sampling cannot be performed.
}

\value{
  A model object inheriting from \code{"gamlss2"} containing posterior samples and
  summary information. The object also contains:
  \itemize{
    \item \code{samples}: Posterior draws for coefficients and smoothing variances.
    \item \code{results}: Model results derived from posterior means.
    \item \code{elapsed}: Elapsed sampling time.
  }
}

\references{
  Umlauf, N., Klein, N., and Zeileis, A. (2018).
  \emph{BAMLSS: Bayesian Additive Models for Location, Scale and Shape (and Beyond)}.
  Journal of Computational and Graphical Statistics, \bold{27}(3), 612--627.
  \doi{10.1080/10618600.2017.1407325}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{bamlss2}}, \code{\link{BS}}
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
## load the abdominal circumference data
data("abdom", package = "gamlss.data")

## specify the model Formula
f <- y ~ s(x) | s(x) | s(x) | s(x)

## estimate model using backfitting
b <- gamlss2(f, data = abdom, family = BCT)

## run MCMC sampling
set.seed(1328)
a <- mcmc(b)

## generate more samples
a <- mcmc(a)

## posterior summary
summary(a)

## plot estimated effects
plot(a)

## plot samples
plot(a, which = "samples")

## compare estimated mean effects
plot(c(b, a))

## predict parameters using samples
pm <- predict(b, FUN = mean)
print(head(pm))
psd <- predict(b, FUN = sd)
print(head(psd))
}

