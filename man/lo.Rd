\name{lo}
\alias{lo}
\alias{special_fit.lo}
\alias{special_predict.lo.fitted}

\title{LOESS Smoothers for GAMLSS}

\description{
  Constructs a local polynomial regression fitting (LOESS) model term for use
  with \code{gamlss2}. The neighborhood size and polynomial degree can be 
  supplied by the user. Neighborhood parameter tuning can be performed 
  automatically using various information criteria during backfitting.
}

\usage{
lo(formula, ...)
}

\arguments{
  \item{formula}{
    A formula specifying the predictors to be included in the local polynomial
    regression (e.g., \code{~ x1} or \code{~ x1 + x2}).
  }
  \item{\dots}{
    Further arguments passed directly to the underlying \code{\link[stats]{loess}} 
    function or its control structures. Typical arguments include \code{span} 
    (the smoothing parameter controlling neighborhood size), \code{degree} 
    (the degree of the local polynomials, usually 1 or 2), and \code{parametric} 
    or \code{drop.square}.
  }
}

\details{
  Function \code{lo()} defines a special model term for \code{gamlss2}. It acts 
  as a wrapper wrapper around the standard \code{\link[stats]{loess}} function 
  integrated cleanly into the iterative backfitting architecture of the location, 
  scale, and shape loops.

  Estimation is based on localized weighted least squares. For a given neighborhood 
  configuration, the local polynomial parameters are updated at every backfitting step 
  using the adjusted working targets \code{z} and weights \code{w}. 

  The effective degrees of freedom (\code{edf}) for the smooth term are extracted 
  analytically using the \code{trace.hat} value returned by the local regression 
  engine.

  To ensure global model identifiability across the additive channels of the 
  distributional parameters, the estimated contribution is centered by subtracting 
  its mean before being passed back to the outer fitting algorithm.
}

\value{
  The constructor \code{lo()} returns an object of classes \code{"special"} 
  and \code{"lo"} containing the covariate matrix, formula definitions, 
  and control arguments required for fitting and prediction.

  After estimation with \code{gamlss2}, the corresponding fitted special term 
  contains the following main components:
  \itemize{
    \item \code{model}: the underlying fitted \code{loess} object.
    \item \code{fitted.values}: centered fitted contribution of the term.
    \item \code{edf}: effective degrees of freedom, taken from \code{trace.hat}.
    \item \code{shift}: the mean value used to center the fitted smooth.
    \item \code{term}: character vector of variable names involved in the smooth.
    \item \code{formula}: the internal updated formula matching the working response.
  }
}

\references{
  Cleveland, W. S., Grosse, E. and Shyu, W. M. (1992).
  \dQuote{Local regression models.}
  Chapter 8 of \emph{Statistical Models in S}, 
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{special_terms}}, \code{\link[stats]{loess}}
}

\examples{\dontrun{## air quality data
air <- subset(airquality, !is.na(Ozone))
air <- air[order(air$Temp), ]

## Fit a distribution where both location (mu) and scale (sigma) 
## depend nonlinearly on Temp via local regression
f <- Ozone ~ lo(~Temp+Wind, span=0.5, degree=2) | lo(Temp, span=0.7)

## estimate model
m <- gamlss2(f, data = air, family = GA)

## model summary
summary(m)

## visualize estimated effects
plot(m)
}}

\keyword{models}
\keyword{smooth}

