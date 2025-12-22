\name{BS} 
\alias{BS}

\title{Bayesian Sampler (BS) for GAMLSS Models}

\description{
  The function \code{BS()} implements an MCMC sampler for generalized additive models for
  location, scale and shape (GAMLSS) as fitted by \code{\link{gamlss2}}. The sampler uses
  a Metropolis-Hastings update for regression coefficients of linear and smooth terms.
  Smoothing variances are updated using univariate slice sampling or conjugate updates,
  depending on the structure of the penalty matrices.
}

\usage{
## Bayesian sampler.
BS(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)
}

\arguments{
  \item{x}{The full model matrix to be used for fitting.}
  \item{y}{The response vector or matrix.}
  \item{specials}{A named list of special model terms, e.g., including design and
    penalty matrices for fitting smooth terms using \code{\link[mgcv]{smooth.construct}}.}
  \item{family}{A family object, see \code{\link{gamlss2.family}}.}
  \item{offsets}{If supplied, a list or data frame of possible model offset.
    See \code{\link{gamlss2}} for details on the offset interface.}
  \item{weights}{If supplied, a numeric vector of prior \code{weights} to be used in the
    fitting process. Should be \code{NULL} or a numeric vector.}
  \item{start}{Starting values for the parameters of the response distribution. See the
    examples for \code{\link{gamlss2}}.}
  \item{xterms}{A named list specifying the linear model terms. Each named list
    element represents one parameter as specified in the family object.}
  \item{sterms}{A named list specifying the special model terms. Each named list
    element represents one parameter as specified in the family object.}
  \item{control}{Further control arguments as specified within the call of \code{\link{gamlss2}}.
    See the details.}
}

\details{
  Function \code{BS()} is typically called via \code{\link{mcmc}} or via wrapper function
  \code{\link{bamlss2}}. It implements a blocked Metropolis-Hastings sampler using
  working responses and working weights. Smooth terms are updated using penalized
  weighted least squares proposals and their smoothing variances are sampled using
  univariate slice sampling.

  The function uses the following control arguments:
  \itemize{
    \item \code{n.iter}: Integer, the total number of MCMC iterations.
    \item \code{burnin}: Integer, the burn-in period.
    \item \code{thin}: Integer, thinning parameter for saved samples.
    \item \code{trace}: Logical, should information be printed while the sampler is running?
    \item \code{flush}: Logical, use \code{\link[utils]{flush.console}} for displaying the
      current output in the console.
    \item \code{nullmodel}: Logical, should an intercept-only null model be evaluated to
      compute a null deviance and deviance reduction?
    \item \code{binning}: Logical, should binning be used for smooth terms if available
      (unique observations for building \eqn{X^{\top}WX} and related quantities)?
  }

  The returned object contains posterior samples for all model components. In addition,
  DIC and related quantities may be computed from saved log-likelihood values if enabled.
}

\value{
  The function returns a fitted model object containing posterior samples and summary
  information. The exact structure is determined by the sampler implementation, but
  typically includes:
  \itemize{
    \item \code{samples}: A matrix of posterior draws for coefficients and smoothing variances.
    \item \code{coefficients}: Posterior means of coefficients.
    \item \code{fitted.values}: Posterior mean fitted values for each distributional parameter.
    \item \code{edf}: Estimated degrees of freedom for smooth terms.
    \item \code{logLik}: Log-likelihood at the posterior mean predictor.
    \item \code{dic}: DIC-related quantities (if computed).
  }
}

\references{
  Umlauf, N., Klein, N., and Zeileis, A. (2018).
  \emph{BAMLSS: Bayesian Additive Models for Location, Scale and Shape (and Beyond)}.
  Journal of Computational and Graphical Statistics, \bold{27}(3), 612--627.
  \doi{10.1080/10618600.2017.1407325}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{RS}}, \code{\link{CG}}, \code{\link{mcmc}}, \code{\link{bamlss2}}
}

\examples{## see ?mcmc
}

