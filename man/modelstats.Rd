\name{Rsq}
\alias{Rsq}
\alias{GAIC}

\title{GAIC and Generalised (Pseudo) R-squared for GAMLSS Models}

\description{
  Functions to compute the GAIC and the generalised R-squared of Nagelkerke (1991) for GAMLSS models.
}

\usage{
## Information criteria.
GAIC(object, ...,
  k = 2, corrected = FALSE)

## R-squared.
Rsq(object, ...,
  type = c("Cox Snell", "Cragg Uhler", "both", "simple"),
  newdata = NULL)
}

\arguments{
  \item{object}{A fitted model object}
  \item{\dots}{Optionally more fitted model objects.}
  \item{k}{Numeric, the \emph{penalty} to be used. The default \code{k = 2} corresponds to
    the classical AIC.}
  \item{corrected}{Logical, whether the corrected AIC should be used?
    Note that it applies only when \code{k = 2}.}
  \item{type}{Which definition of R squared. Can be the \code{"Cox Snell"} or the Nagelkerke, \code{"Cragg Uhler"} or \code{"both"},
    or \code{"simple"}, which computes the R-squared based on the median. In this case also \code{newdata} may be supplied.}
  \item{newdata}{Only for \code{type = "simple"}, the R-squared can be evaluated using \code{newdata}.}
}

\details{
  Note that if \pkg{gamlss} is loaded before \pkg{gamlss2}, the function \code{Rsq()}
  from \pkg{gamlss2} is masked. However, the new \code{Rsq()} method supports fitted
  objects from both \pkg{gamlss} and \pkg{gamlss2}.

  The \code{Rsq()} function uses the definition for R-squared:

  \deqn{R^2=1- \left(\frac{L(0)}{L(\hat{\theta})}\right)^{2/n}}{R-squared =1-(L(0)/L(fitted))^(2/n)}

  where \eqn{L(0)} is the null model (only a constant is  fitted to all parameters) and 
  \eqn{L(\hat{\theta})}{L(fitted)} is the current fitted model. This definition sometimes is
  referred to as the Cox & Snell R-squared. The Nagelkerke /Cragg & Uhler's definition divides
  the above by

  \deqn{1 - L(0)^{2/n}}{1-L(0)^(2/n)}

  Note that \code{GAIC()} function is fully represented by the \code{\link[stats]{AIC}} and
  \code{\link[stats]{BIC}} functions. We only included \code{GAIC()} because of the previous
  package so users who relied on this function would still have the option.
}

\value{
  Numeric vector or data frame, depending on the number of fitted model objects.
}

\references{
  Nagelkerke NJD (1991). \dQuote{A Note on a General Definition of the Coefficient of Determination.}
    \emph{Biometrika}, \bold{78}(3), 691--692. \doi{10.1093/biomet/78.3.691}
}

\seealso{\code{\link{gamlss2}}}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
## load the aids data set
data("aids", package = "gamlss.data")

## estimate negative binomial count models
b1 <- gamlss2(y ~ x + qrt, data = aids, family = NBI)
b2 <- gamlss2(y ~ s(x) + s(qrt, bs = "re"), data = aids, family = NBI)

## compare models
Rsq(b1)
Rsq(b1, type = "both")
Rsq(b1, b2)
GAIC(b1, b2)
AIC(b1, b2)
BIC(b1, b2)
}

\keyword{regression}

