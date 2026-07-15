\name{mcmc}
\alias{mcmc}

\title{Run MCMC Sampling for a Fitted \code{gamlss2} Model}

\description{
  The function \code{mcmc()} runs MCMC sampling for a fitted
  \code{\link{gamlss2}} model. It reuses the fitted model object and merges the
  MCMC output back into the original object so that standard methods such as
  \code{\link{summary}} and \code{\link{plot}} work as expected.
}

\usage{
mcmc(object, n.iter = 1200, burnin = 200, thin = 1)
}

\arguments{
  \item{object}{An object of class \code{"gamlss2"} or \code{"bamlss2"} as
    returned by \code{\link{gamlss2}}.}
  \item{n.iter}{Integer, the total number of MCMC iterations.}
  \item{burnin}{Integer, the burn-in period.}
  \item{thin}{Integer, the thinning parameter for saved samples.}
}

\details{
  The function sets \code{n.iter}, \code{burnin}, and \code{thin} in the model
  control list and obtains starting values from the fitted object, including
  smoothing parameters. The MCMC sampler \code{\link{BS}} is then called and
  the resulting samples and posterior summaries are merged into the original
  model object.

  This approach ensures that the result retains the \code{terms},
  \code{call}, and related components required by extractor and plotting
  methods. If the fitted model object does not contain the design information
  required for MCMC sampling, for example because a reduced control setting was
  used, sampling cannot be performed.
}

\value{
  A model object inheriting from \code{"gamlss2"} containing posterior samples
  and summary information. The object also contains:
  \itemize{
    \item \code{samples}: posterior draws for coefficients and smoothing variances;
    \item \code{results}: model results derived from posterior means;
    \item \code{elapsed}: elapsed sampling time.
  }
}

\references{
  Umlauf N, Klein N, Zeileis A (2018).
  \emph{BAMLSS: Bayesian Additive Models for Location, Scale and Shape (and Beyond)}.
  \emph{Journal of Computational and Graphical Statistics}, \bold{27}(3),
  612--627. \doi{10.1080/10618600.2017.1407325}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{bamlss2}}, \code{\link{BS}}
}

\examples{\dontrun{## specify the model formula
f <- Ozone ~ s(Temp) | . | 1

## estimate model using backfitting
m1 <- gamlss2(f, data = airquality, family = GIG)

## run MCMC sampling
set.seed(1328)
m2 <- mcmc(m1)

## generate more samples
m2 <- mcmc(m2)

## posterior summary
summary(m2)

## plot estimated effects
plot(m2)

## plot samples
plot(m2, which = "samples")
plot(m2, which = "samples", ask = TRUE)
plot(m2, which = "samples", ask = TRUE,
  model = "mu", lag.max = 120, term = "s(Temp)")

## samples can be extracted for further inspection
print(head(m2$samples))

## compare estimated mean effects
plot(c(m1, m2))

## predict parameters using samples
pm <- predict(m2, FUN = mean)
print(head(pm))
psd <- predict(m2, FUN = sd)
print(head(psd))
}}

\keyword{models}
\keyword{regression}
