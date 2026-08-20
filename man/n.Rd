\name{n}
\alias{n}
\alias{special_fit.n}
\alias{special_predict.n.fitted}

\title{Neural Network Model Terms for GAMLSS}

\description{
  Constructs a feed-forward neural network model term for use with
  \code{gamlss2}. The term is fitted using \code{\link[nnet]{nnet}}
  with a linear output unit and is integrated into the iterative
  backfitting algorithm for distributional regression.
}

\usage{
n(formula, ...)
}

\arguments{
  \item{formula}{
    A one-sided formula specifying the predictors entering the neural
    network, for example \code{~ x1 + x2}. Numeric and factor predictors
    can be used.
  }
  \item{\dots}{
    Control arguments for the neural network fit. The currently supported
    arguments are \code{size}, the number of units in the hidden layer;
    \code{maxit}, the maximum number of iterations; \code{decay}, the
    weight-decay parameter; \code{trace}, indicating whether optimization
    progress should be printed; \code{MaxNWts}, the maximum allowable
    number of weights; and \code{scale}, indicating whether numeric
    predictors should be scaled before fitting.

    Defaults are \code{size = 50}, \code{maxit = 1000},
    \code{decay = 0.1}, \code{trace = FALSE},
    \code{MaxNWts = 10000}, and \code{scale = TRUE}.
  }
}

\details{
  Function \code{n()} defines a special neural network model term for
  \code{gamlss2}. At each backfitting update, the current working response
  and weights are passed to \code{\link[nnet]{nnet}} and a feed-forward
  neural network is fitted with \code{linout = TRUE}. Hence, the neural
  network contributes a continuous additive predictor to the corresponding
  distributional parameter.

  The formula supplied to \code{n()} defines the inputs to the network.
  For example, \code{n(~ Temp + Wind)} fits a neural network using
  \code{Temp} and \code{Wind} as inputs.

  By default, each numeric predictor is scaled to the interval \eqn{[0,1]}
  using its observed range before the network is fitted. Factor predictors
  are left unchanged. The same scaling information is stored with the
  fitted term and is applied automatically when predicting from new data.
  Scaling can be disabled with \code{scale = FALSE}.

  The fitted effect is centered by subtracting
  its mean before it is returned to the additive predictor.

  The effective degrees of freedom (\code{edf}) are approximated by the
  number of fitted neural network coefficients returned by
  \code{\link[nnet]{nnet}}.
}

\value{
  The constructor \code{n()} returns an object of classes
  \code{"special"} and \code{"n"} containing the model formula,
  input data, neural network control arguments, and optional scaling
  information required for fitting.

  After estimation with \code{gamlss2}, the corresponding fitted special
  term contains the following main components:
  \itemize{
    \item \code{model}: the underlying fitted \code{nnet} object.
    \item \code{fitted.values}: centered fitted contribution of the
          neural network term.
    \item \code{edf}: effective degrees of freedom, approximated by the
          number of fitted network coefficients.
    \item \code{shift}: the mean value subtracted from the fitted neural
          network contribution.
    \item \code{scalex}: ranges used to scale numeric predictors when
          \code{scale = TRUE}.
    \item \code{transfer}: fitted network weights used to warm-start the
          next backfitting iteration.
  }
}

\references{
  Venables, W. N. and Ripley, B. D. (2002).
  \emph{Modern Applied Statistics with S}, Fourth Edition.
  Springer, New York.
}

\seealso{
  \code{\link{gamlss2}},
  \code{\link{special_terms}},
  \code{\link[nnet]{nnet}}
}

\examples{\dontrun{
## air quality data
air <- subset(airquality, !is.na(Ozone))
air <- air[order(air$Temp), ]

## Neural network effect of Temp and Wind for mu.
## A smaller hidden layer is used here for illustration.
set.seed(1328)
f <- Ozone ~ n(~Temp+Wind,size=50,decay=0.1) | .

## estimate model
m <- gamlss2(f, data = air, family = GA)

## model summary
summary(m)

## visualize estimated effects
plot(m)

## create new data:
## vary temperature while holding wind constant
nd <- with(air, data.frame(
  Temp = seq(min(Temp), max(Temp), length.out = 100),
  Wind = median(Wind)
))

## fitted conditional quantiles
pq <- quantile(
  m,
  newdata = nd,
  probs = c(0.05, 0.5, 0.95)
)

## visualize
plot(Ozone ~ Temp, data = air, pch = 19,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.3),
  xlab = "Temperature (degrees F)",
  ylab = "Ozone (ppb)")
matlines(nd$Temp, pq,
  lwd = 2,
  lty = c(2, 1, 2),
  col = 4)
}}

\keyword{models}
\keyword{neural network}

