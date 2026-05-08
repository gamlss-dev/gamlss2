\name{gamlss2} 
\alias{gamlss2}
\alias{gamlss2.formula}
\alias{gamlss2.list}

\title{Generalized Additive Models for Location Scale and Shape}

\description{
  The function \code{gamlss2()} fits generalized additive models for location,
  scale, and shape (GAMLSS). It estimates one additive predictor for each
  parameter of the chosen response distribution, such as \code{mu},
  \code{sigma}, \code{nu}, or \code{tau}. The number and meaning of these
  parameters are determined by the selected family object, see
  \code{\link{gamlss2.family}}.

  In addition to ordinary linear terms, \code{gamlss2()} supports smooth terms
  from \pkg{mgcv} as well special user-defined model-terms, see \code{\link{specials}}.
}

\usage{
gamlss2(formula, \dots)

\method{gamlss2}{formula}(formula, data, family = NO,
  subset, na.action, weights, offset, start = NULL,
  knots = NULL, control = gamlss2_control(\dots), \dots)

\method{gamlss2}{list}(formula, \dots)
}

\arguments{
  \item{formula}{A GAM-type \code{\link[stats]{formula}} or
    \code{\link[Formula]{Formula}} describing the additive predictors for the
    distribution parameters. All smooth terms from \pkg{mgcv} are supported,
    see also \code{\link[mgcv]{formula.gam}}. For \code{gamlss.list()},
    \code{formula} is a list of \code{\link[stats]{formula}}s.}
  \item{data}{A data frame, list, or environment containing the variables in the model.
    If not found in \code{data}, the variables are taken from \code{environment(formula)},
    typically the environment from which \code{gamlss2()} is called.}
  \item{family}{A \code{gamlss.family} or \code{gamlss2.family} object defining
    the response distribution and the link functions of its parameters.}
  \item{subset}{An optional vector specifying a subset of observations to be used
    in the fitting process.}
  \item{na.action}{\code{NA} processing for setting up the \code{\link[stats]{model.frame}}.}
  \item{weights}{An optional vector of prior \code{weights} to be used in the fitting process.
    Should be \code{NULL} or a numeric vector.}
  \item{offset}{Optional offset terms to be included in the linear predictors during
    fitting. If a single numeric vector is supplied, it is assigned to the first
    parameter of the distribution. If offsets are needed for several parameters,
    provide a data frame or list in the same order as the parameter names in the
    family object.}
  \item{start}{Optional starting values for the estimation algorithm, see
    \code{\link{gamlss2_start}}.}
  \item{knots}{An optional list of user-specified knots, see \code{\link[mgcv]{smoothCon}}.}
  \item{control}{A list of control arguments, see \code{\link{gamlss2_control}}.}
  \item{\dots}{Arguments passed to \code{\link{gamlss2_control}}.}
}

\details{
  The function \code{gamlss2()} is the main entry point for fitting distributional
  regression models in \pkg{gamlss2}.
  \itemize{
    \item The response distribution is specified through a family object, either
      from \pkg{gamlss.dist} or created as a \code{\link{gamlss2.family}} object.
    \item By default, estimation is carried out with the \code{\link{RS}} algorithm.
      The optimizer can be changed through \code{\link{gamlss2_control}}. A family
      object may also provide its own optimizer.
    \item The returned object depends on the chosen optimizer. In the standard case,
      this is an object of class \code{"gamlss2"}, for which summary, plotting,
      prediction, and extractor methods are available.
  }
} 

\value{
  The object returned by the selected optimizer. By default, this is an object
  of class \code{"gamlss2"} produced by \code{\link{RS}}. Standard methods and
  extractor functions are available for this class.
}

\references{
  Rigby RA, Stasinopoulos DM (2005). \dQuote{Generalized Additive Models for Location,
  Scale and Shape (with Discussion).} \emph{Journal of the Royal Statistical Society,
  Series C (Applied Statistics)}, \bold{54}, 507--554.
  \doi{10.1111/j.1467-9876.2005.00510.x}

  Rigby RA, Stasinopoulos DM, Heller GZ, De Bastiani F (2019).
  \emph{Distributions for Modeling Location, Scale, and Shape: Using GAMLSS in R},
  Chapman and Hall/CRC. \doi{10.1201/9780429298547}

  Stasinopoulos DM, Rigby RA (2007). \dQuote{Generalized Additive Models for
  Location Scale and Shape (GAMLSS) in R.} \emph{Journal of Statistical Software},
  \bold{23}(7), 1--46. \doi{10.18637/jss.v023.i07}

  Stasinopoulos DM, Rigby RA, Heller GZ, Voudouris V, De Bastiani F (2017).
  \emph{Flexible Regression and Smoothing: Using GAMLSS in R},
  Chapman and Hall/CRC. \doi{10.1201/b21973}
}

\seealso{
  \code{\link{RS}}, \code{\link{gamlss2_control}}, \code{\link{gamlss2.family}},
  \code{\link{gamlss2_start}}
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
## load the abdominal circumference data
data("abdom", package = "gamlss.data")

## specify a distributional model:
## location and scale vary with x,
## shape parameters remain constant
f <- y ~ s(x) | s(x) | 1 | 1

## estimate model
m1 <- gamlss2(f, data = abdom, family = BCT)

## model summary
summary(m1)

## plot parameter-specific effects
plot(m1, which = "effects")

## plot diagnostics
plot(m1, which = "resid")

## predict selected distribution parameters
par <- predict(m1, model = c("mu", "sigma"))

## fitted centile curves
pq <- quantile(m1, probs = c(0.05, 0.5, 0.95))

## visualize
plot(y ~ x, data = abdom, pch = 19,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
matplot(abdom$x, pq, type = "l", lwd = 2,
  lty = 1, col = 4, add = TRUE)

## count response example using HarzTraffic data
data("HarzTraffic", package = "gamlss2")

## seasonal count model for motorcycles
m2 <- gamlss2(
  bikes ~ s(yday, bs = "cc") | s(yday, bs = "cc"),
  data = HarzTraffic,
  family = NBI
)

## summary
summary(m2)

## effect plots
plot(m2)

## residual diagnostics
plot(m2, which = "resid")

## quantiles over the day of the year
nd <- data.frame(yday = 1:365)
pq <- quantile(m2, newdata = nd, probs = c(0.05, 0.5, 0.95))

## visualize fitted quantiles
plot(bikes ~ yday, data = HarzTraffic, pch = 19,
  col = gray(0.1, alpha = 0.3))
matplot(nd$yday, pq, type = "l", lty = c(2, 1, 2),
  lwd = 2, col = 4, add = TRUE)
}

\keyword{models}
\keyword{regression}
