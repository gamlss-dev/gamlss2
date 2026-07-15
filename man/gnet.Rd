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
  \item{formula}{A formula specifying the covariates that should be estimated using
    the \pkg{glmnet} implementation.}
  \item{\dots}{Control arguments passed to \code{\link[glmnet]{glmnet}} and the
    information-criterion-based shrinkage selection.}
}

\details{
  \code{gnet()} constructs a \code{\link[stats]{model.matrix}} representation
  of the supplied formula and delegates estimation to \code{\link[glmnet]{glmnet}}
  inside the backfitting algorithm. The optimal penalty parameter is selected
  using an information criterion, with \code{criterion = "bic"} used by default.
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
sp <- specials(m, model = "mu", term = "gnet(", element = "model")
print(class(sp))
plot(sp)
}}
