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
