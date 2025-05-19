\name{find_family}
\alias{find_family}
\alias{fit_family}
\alias{available_families}
\alias{find_gamlss2}

\title{Find and Fit GAMLSS Families}

\description{
  These functions provide useful infrastructures for finding suitable GAMLSS
  families for a response variable.
}

\usage{
## List of available families from gamlss.dist package.
available_families(type = c("continuous", "discrete"), families = NULL, ...)

## Find suitable response distribution.
find_family(y, families = NULL, k = 2, verbose = TRUE, ...)

## Fit distribution parameters.
fit_family(y, family = NO, plot = TRUE, ...)

## Find gamlss2 model.
find_gamlss2(formula, families = NULL, k = 2,
  select = FALSE, verbose = TRUE, ...)
}

\arguments{
  \item{type}{Character, is the reponse continuous or discrete?}
  \item{families}{Character, the names of the family objects of the \pkg{gamlss.dist}
    package that should be returned.}
  \item{y}{The response vector or matrix.}
  \item{k}{Numeric, the penalty factor that should be used for the \code{\link{GAIC}}.}
  \item{verbose}{Logical, should runtime information be printed?}
  \item{family}{A famnily object that should be used for estimation,
    see also \code{\link{gamlss2.family}}.}
  \item{plot}{Logical, should a plot of the fitted density be provided?}
  \item{formula}{A model formula, see \code{\link{gamlss2}}.}
  \item{select}{Logical, if set to \code{select = TRUE}, model term selection is
    enforced using model fitting function \code{\link{select_gamlss2}}.}
  \item{\dots}{For function \code{available_families()}, arguments passed to the family objects,
    e.g., for setting link functions \code{mu.link = "log"}.
    Further arguments to be passed to \code{\link{gamlss2}} when using
    \code{find_family()}, or arguments \code{legend = TRUE/FALSE}, \code{pos = "topright"} (see also
    function \code{\link[graphics]{legend}}),
    \code{main}, \code{xlab} and \code{ylab} when argument \code{plot = TRUE} using
    function \code{fit_family()}. For function \code{find_gamlss()}, arguments
    are passed to \code{available_families()} and \code{\link{gamlss2}}.}
}

\details{
  The function \code{find_family()} employs \code{\link{gamlss2}} to estimate intercept-only
  models for each specified family object in the \code{families} argument. Note that model
  estimation occurs within a \code{\link[base]{try}} block with warnings suppressed.
  Additionally, the function calculates the \code{\link{GAIC}} for each family whenever
  feasible and returns the sorted values in descending order.

  Function \code{fit_family()} fits a single intercept-only model using the
  specified family and creates a plot of the fitted density.
}

\value{
  Function \code{find_family()} returns a vector of \code{\link{GAIC}} values for the
  different fitted families. Function \code{fit_family()} returns the fitted intercept-only
  model.
}

\seealso{
  \code{\link{gamlss2}}.
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
\dontrun{## load data
data("rent", package = "gamlss.data")

## find a suitable response to the response
ic <- find_family(rent$R)
print(ic)

## fit parameters using the BCCG family
fit_family(rent$R, family = BCCG)

## count data
data("polio", package = "gamlss.data")

## search best count model
ic <- find_family(polio, k = 0,
  families = available_families(type = "discrete"))
print(ic)

## fit parameters using the ZASICHEL family
fit_family(polio, family = ZASICHEL)

## search complete model
## b <- find_gamlss2(R ~ s(Fl) + s(A) | . | . | ., data = rent,
##   select = TRUE, mu.link = "log")
}
}

\keyword{regression}
\keyword{models}
\keyword{distribution}

