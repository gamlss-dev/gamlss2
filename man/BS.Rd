\name{BS}
\alias{BS}

\title{Bayesian Sampler for GAMLSS Models}

\description{
  The function \code{BS()} implements an MCMC sampler for generalized additive
  models for location, scale, and shape (GAMLSS) as fitted by
  \code{\link{gamlss2}}. The sampler uses Metropolis-Hastings updates for
  regression coefficients of linear and smooth terms.
}

\usage{
BS(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)
}

\arguments{
  \item{x}{The full model matrix to be used for fitting.}
  \item{y}{The response vector or matrix.}
  \item{specials}{A named list of special model terms, for example including
    design and penalty matrices for fitting smooth terms using
    \code{\link[mgcv]{smooth.construct}}.}
  \item{family}{A family object, see \code{\link{gamlss2.family}}.}
  \item{offsets}{If supplied, a list or data frame of model offsets. See
    \code{\link{gamlss2}} for details on the offset interface.}
  \item{weights}{If supplied, a numeric vector of prior weights used in the
    fitting process.}
  \item{start}{Starting values for the parameters of the response distribution.
    See the examples for \code{\link{gamlss2}}.}
  \item{xterms}{A named list specifying the linear model terms. Each named list
    element represents one parameter as specified in the family object.}
  \item{sterms}{A named list specifying the special model terms. Each named list
    element represents one parameter as specified in the family object.}
  \item{control}{Further control arguments as specified in the call to
    \code{\link{gamlss2}}.}
}

\details{
  Function \code{BS()} is typically called via \code{\link{mcmc}} or via the
  wrapper function \code{\link{bamlss2}}. It implements a blocked
  Metropolis-Hastings sampler using working responses and working weights.
  Smooth terms are updated using penalized weighted least squares proposals and
  their smoothing variances are sampled using univariate slice sampling.

  The function uses the following control arguments:
  \itemize{
    \item \code{n.iter}: integer, the total number of MCMC iterations;
    \item \code{burnin}: integer, the burn-in period;
    \item \code{thinning}: integer, the thinning parameter for saved samples;
    \item \code{trace}: logical, should information be printed while the sampler is running;
    \item \code{flush}: logical, whether to use \code{\link[utils]{flush.console}}
      to display current output in the console;
    \item \code{nullmodel}: logical, should an intercept-only null model be
      evaluated to compute a null deviance and deviance reduction?;
    \item \code{binning}: logical, should binning be used for smooth terms if available?
  }
}

\value{
  A fitted model object containing posterior samples and summary information.
  Typically it includes:
  \itemize{
    \item \code{samples}: posterior draws for coefficients and smoothing variances;
    \item \code{coefficients}: posterior means of coefficients;
    \item \code{fitted.values}: posterior mean fitted values for each distributional parameter;
    \item \code{edf}: estimated degrees of freedom for smooth terms;
    \item \code{logLik}: log-likelihood at the posterior mean predictor.
  }
}

\references{
  Umlauf N, Klein N, Zeileis A (2018).
  \emph{BAMLSS: Bayesian Additive Models for Location, Scale and Shape (and Beyond)}.
  \emph{Journal of Computational and Graphical Statistics}, \bold{27}(3),
  612--627. \doi{10.1080/10618600.2017.1407325}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{RS}}, \code{\link{CG}},
  \code{\link{mcmc}}, \code{\link{bamlss2}}
}

\examples{
## see ?mcmc
}

\keyword{models}
\keyword{regression}
