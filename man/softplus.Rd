\name{softplus}
\alias{softplus}

\title{
  Softplus Link Object
}

\description{
  Implementation of the softplus link function.
}

\usage{
softplus(a = 1)
}

\arguments{
  \item{a}{Extra parameter of the generalized softplus function}
}

\value{
  An object of class \code{"link-glm"}.
}

\seealso{
  \code{\link{gamlss2}}.
}

\examples{
## softplus function for different values of a.
plot(softplus(a = 1)$linkinv, -4, 4, ylab = "softplus(x)")
plot(softplus(a = 10)$linkinv, -4, 4, col = 4, add = TRUE)
legend("topleft", c("a = 1", "a = 10"), lwd = 1, col = c(1, 4))
}

\keyword{regression}
\keyword{models}
\keyword{distribution}

