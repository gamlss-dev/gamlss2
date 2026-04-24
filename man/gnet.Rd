\name{gnet}
\alias{gnet}

\title{Lasso with \pkg{glmnet}}

\description{
  Compatibility constructor for estimating Lasso model terms using
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
  \code{gnet()} is a compatibility wrapper around the current
  \code{\link{lasso}}/\code{"glmnet"} special-term path. It constructs a
  \code{\link[stats]{model.matrix}} representation and delegates estimation
  to \code{\link[glmnet]{glmnet}} inside the backfitting algorithm.
}

\value{
  A special model-term object used internally by \code{\link{gamlss2}}.
}

\seealso{
  \code{\link{lasso}}, \code{\link{la}}, \code{\link{gamlss2}}, \code{\link{specials}}.
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
\dontrun{data("rent", package = "gamlss.data")

rent$Flc <- cut(rent$Fl, breaks = seq(20, 160, by = 10),
  include.lowest = TRUE)
rent$Ac <- cut(rent$A, breaks = seq(1890, 1990, by = 10),
  include.lowest = TRUE)

f <- R ~ gnet(~Flc + Ac + loc) |
  gnet(~Flc + Ac + loc) |
  gnet(~Flc + Ac + loc) |
  gnet(~Flc + Ac + loc)

m <- gamlss2(f, data = rent, family = BCT)

summary(m)
}}
