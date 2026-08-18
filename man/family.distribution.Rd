\name{family.distribution}
\alias{family.distribution}
\alias{family.GAMLSS}
\alias{family.gamlss.family}
\alias{distributions3_family}

\title{Generate a 'gamlss2' Family from a 'distributions3' Object}

\description{
Set up a [gamlss2.family] object based on an object
inheriting from `distribution` (as provided by the \pkg{distributions3}
package).
}

\usage{
\method{family}{distribution}(object, links, score = TRUE, hessian = FALSE, update = FALSE, ...)

\method{family}{GAMLSS}(object, ...)

\method{family}{gamlss.family}(object, score = TRUE, hessian = FALSE, update = FALSE, ...)

distributions3_family(
  distribution,
  links,
  score = TRUE,
  hessian = FALSE,
  update = FALSE,
  create_distribution = "structure",
  type = NULL,
  valid.response = NULL,
  initialize = NULL,
  ...
)
}

\arguments{
\item{object}{an object specifying a distribution, typically created by
\code{\link[distributions3]{Normal}} or \code{\link[distributions3]{Binomial}}
or \code{\link[gamlss.dist]{GAMLSS}}. Additionally, a
\code{\link[gamlss.dist]{gamlss.family}} from the \pkg{gamlss.dist} package
is also supported.}

\item{links}{named vector or list. Specification of one link per parameter to
be modeled, either via a character string to be passed to \code{\link{make.link2}} or
a \code{link-gamlss2} object (inheriting from \code{link-glm}).}

\item{score}{logical. Should the \code{\link[distributions3]{score}} method for
the distribution object be used to compute the scores? If set to \code{FALSE} they are
approximated numerically.}

\item{hessian}{logical. Should the \code{\link[distributions3]{hessian}} method for
the distribution object be used to compute the Hessians? If set to \code{FALSE} they are
approximated numerically.}

\item{update}{logical. Should the \code{\link[distributions3]{score}} and
\code{\link[distributions3]{hessian}} methods for the distribution object be used
to compute the updates of the linear predictor and the corresponding weights?
If set to \code{FALSE} they are approximated numerically.}

\item{...}{arguments passed from the methods to the workhorse function
\code{distributions3_family}.}

\item{distribution}{character. Name of the class and distribution object creator,
e.g., \code{"Normal"} or \code{"Binomial"}.}

\item{create_distribution}{either a function or one of the character strings
\code{"structure"} or \code{"do.call"}. This specifies how to create a distribution object
from a (fitted) parameter list (typically from [gamlss2]). Either a dedicated
function can be specified or such a function is set up internally. If
\code{create_distribution = "structure"} then only the class attribute is added in
a \code{structure()} call. If set to \code{"do.call"} then the distribution object is
created via \code{do.call()} to the class creator (such as \code{\link[distributions3]{Normal}}
etc.). The latter is safer but requires much more overhead.}

\item{type}{character. Specification of the type of distribution (\code{"continuous"}
vs. \code{"discrete"}).}

\item{valid.response}{function. A \code{function(x)} returning a logical indicating whether
\code{x} is a valid response for the distribution.}

\item{initialize}{list of functions. A named list with one \code{function(y, ...)}
per parameter to be modeled, returning an initialization for this parameter.}
}

\details{
The workhorse function \code{distributions3_family} sets up a \code{\link{gamlss2.family}}
list by leveraging the generic functions and methods from the \pkg{distributions3}
package. The main idea is that \code{\link{gamlss2}} in each iteration and in the final
model fit provides a list of fitted parameters. These can then be mapped to the corresponding
distribution class and all the corresponding methods (such as \code{\link[distributions3]{pdf}}
for computing the probability density function) from \pkg{distributions3} can be leveraged
to compute the building blocks for fitting the model (e.g., \code{pdf}, \code{score},
\code{hessian}, etc.) and for assessing the final model fit.

The \code{family} methods all call \code{distributions3_family} internally but help
to automatically infer and set up some of the elements of the family.
}

\value{
A list inheriting from \code{\link{gamlss2.family}}.
}

\examples{
library("distributions3")
f <- family(Poisson(), links = c(lambda = "log"))
str(f)
}
