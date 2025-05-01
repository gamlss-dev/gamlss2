\name{GDF}

\alias{GDF}
\alias{format.GDF}
\alias{print.GDF}
\alias{mean.GDF}
\alias{variance.GDF}
\alias{skewness.GDF}
\alias{kurtosis.GDF}
\alias{pdf.GDF}
\alias{log_pdf.GDF}
\alias{cdf.GDF}
\alias{quantile.GDF}
\alias{random.GDF}
\alias{support.GDF}
\alias{is_discrete.GDF}
\alias{is_continuous.GDF}

\title{Create a GDF Distribution}

\description{
The GDF (gamlss2 Distribution Family) is a unified class with corresponding methods
that represent all distributional families supported by the \pkg{gamlss2} package. It 
enables seamless integration with the \pkg{distributions3} workflow and provides a 
consistent interface for model fitting and distributional computations.
}

\usage{
GDF(family, parameters)
}

\arguments{
  \item{family}{character. Name of a \code{\link{gamlss2.family}} or a family provided
    by the \pkg{gamlss.dist} package, e.g,
    \code{\link{NO}} or \code{\link{BI}} for the normal or
    binomial distribution, respectively.}
  \item{parameters}{numeric, matrix, list or data frame, see the examples.}
}

\details{
The S3 class \code{GDF} is a slightly more general implementation of the S3 class
\code{\link[gamlss.dist]{GAMLSS}}. For details please see the documentation of \code{\link[gamlss.dist]{GAMLSS}}
}

\value{
A \code{GDF} object, inheriting from \code{distribution}.
}

\references{
Zeileis A, Lang MN, Hayes A (2022).
\dQuote{distributions3: From Basic Probability to Probabilistic Regression.}
Presented at \emph{useR! 2022 - The R User Conference}.
Slides, video, vignette, code at \url{https://www.zeileis.org/news/user2022/}.
}

\seealso{
\code{\link{gamlss2.family}}
}

\examples{
\dontshow{ if(!requireNamespace("distributions3")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }
}
## package and random seed
library("distributions3")
set.seed(6020)

## one normal distribution
X <- GDF("NO", c(mu = 1, sigma = 2))
X

## two normal distributions
X <- GDF("NO", cbind(c(1, 1.5), c(0.6, 1.2)))
X

## three Weibull distributions
X <- GDF("WEI", list(mu = c(1, 1, 2), sigma = c(1, 2, 2)))
X

## see ?gamlss.dist::GAMLSS for the remainder of this example
}

