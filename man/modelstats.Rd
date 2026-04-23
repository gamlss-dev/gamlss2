\name{Rsq}
\alias{Rsq}
\alias{GAIC}

\title{GAIC and Generalised (Pseudo) R-squared for GAMLSS Models}

\description{
  Functions to compute the GAIC and the generalized R-squared of Nagelkerke (1991)
  for GAMLSS models.
}

\usage{
## information criteria
GAIC(object, ...,
  k = 2, corrected = FALSE)

## r-squared
Rsq(object, ...,
  type = c("Cox Snell", "Cragg Uhler", "both", "simple"),
  newdata = NULL)
}

\arguments{
  \item{object}{A fitted model object.}
  \item{\dots}{Optionally more fitted model objects.}
  \item{k}{Numeric, the \emph{penalty} to be used. The default \code{k = 2} corresponds to
    the classical AIC.}
  \item{corrected}{Logical, whether the corrected AIC should be used?
    Note that it applies only when \code{k = 2}.}
  \item{type}{Which definition of R-squared should be used?
    Possible values are \code{"Cox Snell"}, \code{"Cragg Uhler"},
    \code{"both"}, and \code{"simple"}. The option \code{"simple"}
    computes an R-squared based on the median; in this case \code{newdata}
    may also be supplied.}
  \item{newdata}{For \code{type = "simple"}, the R-squared can also be evaluated on
    \code{newdata}.}
}

\details{
  The \code{Rsq()} function uses the following definition of R-squared:

  \deqn{R^2=1- \left(\frac{L(0)}{L(\hat{\theta})}\right)^{2/n}}{R-squared =1-(L(0)/L(fitted))^(2/n)}

  where \eqn{L(0)} is the null model (only a constant is fitted to all parameters) and
  \eqn{L(\hat{\theta})}{L(fitted)} is the fitted model under consideration. This definition is
  sometimes referred to as the Cox and Snell R-squared. The Nagelkerke or Cragg and Uhler
  definition divides the above by

  \deqn{1 - L(0)^{2/n}}{1-L(0)^(2/n)}
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
m1 <- gamlss2(y ~ x + qrt, data = aids, family = NBI)
m2 <- gamlss2(y ~ s(x) + s(qrt, bs = "re"), data = aids, family = NBI)

## compare models
Rsq(m1)
Rsq(m1, type = "both")
Rsq(m1, m2)
GAIC(m1, m2)
AIC(m1, m2)
BIC(m1, m2)

## plot estimated effects
plot(m2)
}

\keyword{regression}
