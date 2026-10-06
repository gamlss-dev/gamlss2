\name{gnet}
\alias{gnet}

\title{Model Terms with \pkg{glmnet}}

\description{
  Constructor for estimating penalized model terms using
  the \pkg{glmnet} package.
}

\usage{
gnet(formula, ...)
}

\arguments{
  \item{formula}{A formula specifying the covariates for the penalized term.}
  \item{\dots}{\code{criterion} selects the information criterion:
    \code{"bic"} (default), \code{"aic"}, \code{"aicc"}, \code{"gaic"}, or
    \code{"gcv"}. Other arguments are passed to \code{\link[glmnet]{glmnet}},
    including \code{alpha} (elastic-net mixing parameter), \code{lambda}
    (penalty sequence), \code{nlambda} (number of penalties), and
    \code{standardize} (whether to standardize covariates).}
}

\details{
  \code{gnet()} constructs a \code{\link[stats]{model.matrix}} from the formula,
  removes its intercept column, and fits the term using
  \code{\link[glmnet]{glmnet}} within backfitting. The penalty \code{lambda}
  minimizes the selected criterion. For \code{criterion = "gaic"}, the
  multiplier \code{K} is set in \code{\link{gamlss2_control}} (default: 2).
}

\value{
  A special model-term object used internally by \code{\link{gamlss2}}.
}

\seealso{
  \code{\link{la}}, \code{\link{gamlss2}}, \code{\link{specials}}.
}

\examples{\dontrun{## transform to categorical
air <- subset(airquality, !is.na(Ozone))

air$Tc <- cut(air$Temp, breaks = seq(50, 100, by = 10),
  include.lowest = TRUE)
air$Wc <- cut(air$Wind, breaks = seq(2, 22, by = 4),
  include.lowest = TRUE)

## model formula, for each parameter we use
## the same gnet() model term
f <- Ozone ~ gnet(~Tc + Wc) | .

## estimate model
m <- gamlss2(f, data = air, family = NO)

## summary with edf
summary(m)

## plot coefficient paths
sp <- specials(m, parameter = "mu", terms = "gnet(~Tc+Wc)", elements = "model")
print(class(sp))
plot(sp)
}}

