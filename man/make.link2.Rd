\name{make.link2}
\alias{make.link2}

\title{Create a Link for Families}

\description{
  This function is used with the \code{family} functions in \code{gamlss2()}.
  Given the name of a link, it returns a link function, an inverse
  link function, the derivative \eqn{d\mu/d\eta}, and a function for domain
  checking. Note that \code{make.link2()} is slightly more flexible and also allows
  functions as arguments.
}

\usage{
make.link2(link)
}

\arguments{
  \item{link}{
    A character string, see function \code{\link[stats]{make.link}}, or function.
  }
}

\value{
  A list containing the following components:
  \item{linkfun}{Link function \code{function(mu)}.}
  \item{linkinv}{Inverse link function \code{function(eta)}.}
  \item{mu.eta}{Derivative \code{function(eta)}: \eqn{d\mu/d\eta}.}
  \item{valideta}{Function \code{function(eta)} that returns \code{TRUE} if \code{eta} is in the domain of \code{linkinv}.}
  \item{name}{A character string representing the name of the link function.}
}

\seealso{
  \code{\link[stats]{make.link}}, \code{\link{gamlss2}}, \code{\link{gamlss2.family}}.
}

\examples{
## Character specification.
utils::str(make.link2("logit"))

## Functions.
utils::str(make.link2(softplus))
}

