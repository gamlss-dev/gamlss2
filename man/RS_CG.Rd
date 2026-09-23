\name{RS} 
\alias{RS}
\alias{CG}

\title{Rigby and Stasinopoulos (RS) & Cole and Green (CG) Algorithm}

\description{
  The function \code{RS()} implements the algorithm of Rigby and Stasinopoulos,
  the function \code{CG()} implements the algorithm of Cole and Green for estimating
  a GAMLSS with \code{\link{gamlss2}}.
}

\usage{
## rigby and stasinopoulos algorithm
RS(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)

## cole and green algorithm
CG(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)
}

\arguments{
  \item{x}{The full model matrix to be used for fitting.}
  \item{y}{The response vector or matrix.}
  \item{specials}{A named list of special model terms, e.g., including design and
    penalty matrices for fitting smooth terms using \code{\link[mgcv]{smooth.construct}}.}
  \item{family}{A family object, see \code{\link{gamlss2.family}}.}
  \item{offsets}{If supplied, a list or data frame of possible model offset.}
  \item{weights}{If supplied, a numeric vector of weights.}
  \item{start}{Starting values, either for the parameters of the response distribution or,
    if specified as a named list in which each element of length one is named with
    \code{"(Intercept)"}, the respective intercepts are initialized. If starting values are
    specified as a named list, data frame or matrix, where each element/column is a
    vector with the same length as the number of observations in the data,
    the respective predictors are initialized with these. See the examples for
    \code{\link{gamlss2}}.}
  \item{xterms}{A named list specifying the linear model terms. Each named list
    element represents one parameter as specified in the family object.}
  \item{sterms}{A named list specifying the special model terms. Each named list
    element represents one parameter as specified in the family object.}
  \item{control}{Further control arguments as specified within the call of \code{\link{gamlss2}}.
    See the details.}
}

\details{
  Functions \code{RS()} and \code{CG()} are called within \code{\link{gamlss2}}. Both functions
  implement a backfitting algorithm for estimating GAMLSS. For algorithm details see
  Rigby and Stasinopoulos (2005).

  The functions use the following control arguments:
  \itemize{
    \item \code{eps}: Numeric vector of length 1 to 3, the stopping criterion
      for the global, per-parameter, and CG correction loops. The CG tolerance
      defaults to the per-parameter tolerance.
    \item \code{maxit}: Integer vector of length 1 to 3, the maximum number of global
      parameter sweeps and per-parameter IRLS/backfitting iterations. Default is
      \code{maxit = c(20, 50)}. For the CG algorithm, an optional third element
      limits the internal correction sweeps and defaults to the second element.
      Plain RS always performs
      one parameter sweep before checking global convergence; a third element is
      therefore ignored by RS.
    \item \code{step}: Numeric, the step length control parameter. Default is
      \code{step = 1}. After the first update and with fixed smoothing
      parameters, the step is applied to a complete parameter update in RS and
      to a complete correction sweep in CG. A smoothing-parameter change and
      its coefficients are instead accepted or rejected as one state. If
      \code{step} is set smaller than 1, it might be appropriate to lower the
      stopping criterion \code{eps}, too.
    \item \code{autostep}: Logical, should an update that decreases its
      objective be halved automatically? Default is \code{TRUE}. The objective
      is the likelihood for unpenalized updates and the penalized likelihood
      when smoothing parameters are fixed. Smoothing-parameter changes use the
      specified smoothing criterion, the default local REML update, or the
      term-specific local ML update.
    \item \code{sigma.tol}: Numeric between zero and one. If a model contains
      a smooth for \code{sigma}, warn when its smallest fitted value divided by
      its median is less than \code{sigma.tol}. The default is \code{0.001};
      use zero or \code{FALSE} to disable the diagnostic. This diagnostic does
      not constrain the fitted values.
    \item \code{CG}: Integer, the outer iteration in which to start the CG
      correction. Thus, \code{CG = 5} uses RS for the first four outer
      iterations and CG in the fifth. \code{CG = 0} starts with CG immediately;
      the default is \code{CG = Inf}.
    \item \code{trace}: Logical, should information be printed while the algorithm is running?
    \item \code{rs.cache}: Logical, reuse unchanged generated family evaluations
      and weighted matrix calculations within a fit? Default is \code{TRUE}.
      Reuse requires identical inputs and preserves the original arithmetic.
      User replacements of family callbacks are not memoized. Linear QR reuse
      is restricted to full-rank fits with strictly positive weights.
      Set \code{FALSE} to disable these caches, for example for comparisons.
    \item \code{flush}: Logical, use \code{\link[utils]{flush.console}} for displaying the
      current output in the console.
    \item \code{ridge}: Logical, should automatic ridge penalization be applied
      only to linear effects, without penalizing the intercept? For each parameter of
      the distribution the optimum ridge penalty is estimated using an information criterion.
      Possible options are \code{criterion = c("aic", "aicc", "bic", "gaic", "gcv")}. The
      default is \code{criterion = "gaic"} and argument \code{K = 2}, which can
      be set in \code{\link{gamlss2_control}}.
  }

  To facilitate the development of new algorithms for \pkg{gamlss2}, users can exchange
  them using the \code{optimizer} argument in \code{\link{gamlss2_control}}. Users
  developing new model fitting functions are advised to use these functions as templates
  and pass them to \code{\link{gamlss2_control}}. Alternatively, users can replace the
  optimizer function by adding a named list element, \code{"optimizer"}, to the family object.
  For instructions on setting up new families in \pkg{gamlss2}, see \code{\link{gamlss2.family}}.
} 

\value{
  Functions \code{RS()} and \code{CG()} return a named list of class \code{"gamlss2"}
  containing the following objects:
  \item{fitted.values}{A data frame of the fitted values of the predictors
    of the selected distribution.}
  \item{fitted.specials}{A named list, one element for each parameter of the distribution,
    containing the fitted model object information of special model terms.}
  \item{fitted.linear}{A named list, one element for each parameter of the distribution,
    containing the information on fitted linear effects.}
  \item{coefficients}{A named list, one element for each parameter of the distribution,
    containing the estimated parameters of the linear effects.}
  \item{elapsed}{The elapsed runtime of the algorithm.}
  \item{iterations}{How many iterations the algorithm performed.}
  \item{converged}{Logical, did the global iteration satisfy the stopping criterion?}
  \item{logLik}{The final value of the log-likelihood of the model.}
  \item{control}{All effective control arguments used by the optimizer.}
}

\references{
  Rigby RA, Stasinopoulos DM (2005). \dQuote{Generalized Additive Models for Location,
  Scale and Shape (with Discussion).} \emph{Journal of the Royal Statistical Society,
  Series C (Applied Statistics)}, \bold{54}, 507--554.
  \doi{10.1111/j.1467-9876.2005.00510.x}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2_control}}, \code{\link{gamlss2.family}}
}

\examples{\dontrun{## count response example using HarzTraffic data
data("HarzTraffic", package = "gamlss2")

## specify the model formula
f <- bikes ~ s(yday, bs = "cc") + s(temp) + s(rain) | .

## estimate model using RS (default)
m1 <- gamlss2(f, data = HarzTraffic, family = NBI, optimizer = RS)

## now with CG
m2 <- gamlss2(f, data = HarzTraffic, family = NBI, optimizer = CG)

## first 2 RS iterations and afterwards switch to CG
m3 <- gamlss2(f, data = HarzTraffic, family = NBI, CG = 3)
}}

\keyword{models}
\keyword{regression}
