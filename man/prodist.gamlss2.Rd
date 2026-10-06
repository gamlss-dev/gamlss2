\name{prodist.gamlss2}

\alias{prodist.gamlss2}

\title{Extracting Fitted or Predicted Probability Distributions from gamlss2 Models}

\description{
  Methods for \pkg{gamlss2} model objects for extracting fitted (in-sample) or
  predicted (out-of-sample) probability distributions as \pkg{distributions3}
  objects.
}

\usage{
\method{prodist}{gamlss2}(object, ...)
}
\arguments{
  \item{object}{A model object of class \code{\link{gamlss2}}.}
  \item{...}{Arguments passed on to \code{\link{predict.gamlss2}}, e.g., \code{newdata}.}
}

\details{
  For ordinary predictions, \code{prodist(object, newdata = ...)} returns the
  same distributions as the default \code{predict(object, newdata = ...)}.
  Use \code{\link[distributions3]{prodist}} when code needs a common interface
  for extracting probability distributions from different supported model
  classes. Omit \code{newdata} to obtain fitted distributions.

  The returned objects support \pkg{distributions3} methods such as
  \code{pdf()}, \code{cdf()}, \code{quantile()}, and \code{random()}.
  They describe response variation conditional on the fitted parameters;
  they do not include parameter-estimation uncertainty.
}

\value{
  A \pkg{distributions3} object inheriting from \code{distribution}.
  Its concrete class depends on the family: a \pkg{distributions3} family
  retains its distribution class (e.g., \code{Normal}); families without a
  distribution constructor use \code{\link{GAMLSS2}}.
}

\references{
Zeileis A, Lang MN, Hayes A (2022).
\dQuote{distributions3: From Basic Probability to Probabilistic Regression.}
Presented at \emph{useR! 2022 - The R User Conference}.
Slides, video, vignette, code at \url{https://www.zeileis.org/news/user2022/}.
}

\seealso{
\code{\link{GAMLSS2}}, \code{\link{family.distribution}},
\code{\link{predict.gamlss2}}
}

\examples{\dontrun{## fit a Gamma model for stopping distance given speed
library("distributions3")
data("cars", package = "datasets")
m <- gamlss2(dist ~ s(speed) | s(speed),
  data = cars, family = GA)

## prodist() and default predict() return the same distributions
nd <- data.frame(speed = c(10, 20, 30))
d <- prodist(m, newdata = nd)
all.equal(d, predict(m, newdata = nd))
print(d)

## evaluate the predicted distributions
quantile(d, c(0.05, 0.5, 0.95), elementwise = FALSE)
mean(d)
variance(d)
pdf(d, c(10, 50, 100))
cdf(d, c(10, 50, 100))
random(d, 5)

## fitted distributions reproduce the model log-likelihood
sum(log_pdf(prodist(m), cars$dist))
logLik(m)
}}

\keyword{distribution}

