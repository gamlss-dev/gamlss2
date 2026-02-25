\name{lasso}
\alias{lasso}

\title{Lasso with \pkg{glmnet}}

\description{
  Constructor function for estimating Lasso model terms using
  the \pkg{glmnet} package.
}

\usage{
## Model term constructor function.
lasso(formula, ...)
}

\arguments{
  \item{formula}{A formula specifying the covariates that should be estimated using
    the Lasso implementation in \pkg{glmnet}.}
  \item{\dots}{Control arguments to be passed to function \code{\link[glmnet]{glmnet}}.}
}

\details{
  The \code{formula} is used to generate a \code{\link[stats]{model.matrix}},
  which is then used for estimation. Note that the data is not scaled, so the user
  must scale the covariates manually. The function \code{\link[glmnet]{glmnet}} is
  employed within the backfitting algorithm implemented in \code{\link{RS}}.
  The optimal Lasso shrinkage parameter is selected based on an information
  criterion. Available options for the criterion are
  \code{criterion = c("gcv", "aic", "gaic", "aicc", "bic")}.
}

\value{
  The \code{lasso()} function is used internally within \code{\link{gamlss2}}
  and provides the necessary details for estimating Lasso-type model terms using
  function \code{\link[glmnet]{glmnet}}.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{specials}}.
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
\dontrun{## load data
data("rent", package = "gamlss.data")

## transform numeric to factor variables
rent$Flc <- cut(rent$Fl, breaks = seq(20, 160, by = 10),
  include.lowest = TRUE)
rent$Ac <- cut(rent$A, breaks = seq(1890, 1990, by = 10),
  include.lowest = TRUE)

## set up the model formula for a BCT model
f <- R ~ lasso(~Flc + Ac + loc) |
  lasso(~Flc + Ac + loc) |
  lasso(~Flc + Ac + loc) |
  lasso(~Flc + Ac + loc)

## estimation
b <- gamlss2(f, data = rent, family = BCT)

## summary, shows the estimated degrees of freedom
## for each model term
summary(b)

## extract fitted special lasso model term
st <- specials(b, model = "mu")

## plot coefficient paths
plot(st$model)

## same for sigma
plot(specials(b, model = "sigma", elements = "model"))
}}

\keyword{regression}
\keyword{models}

