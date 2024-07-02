\name{gamlss2-package}
\alias{gamlss2-package}

\docType{package}

\title{GAMLSS Modeling with Advanced Flexible Infrastructures}

\description{
  \packageDescription{gamlss2}
}

\details{
  The primary purpose of this package is to facilitate the creation of advanced infrastructures
  designed to enhance the Generalized Additive Models for Location Scale and Shape (GAMLSS,
  Rigby and Stasinopoulos 2005) modeling framework. Notably, the \pkg{gamlss2} package represents a
  significant overhaul of its predecessor, \pkg{gamlss}, with a key emphasis on improving estimation
  speed and incorporating more adaptable infrastructures. These enhancements enable the seamless
  integration of various algorithms into GAMLSS, including gradient boosting, Bayesian estimation,
  regression trees, and forests, fostering a more versatile and powerful modeling environment.

  Moreover, the package expands its compatibility by supporting all model terms from the base
  R \pkg{mgcv} package. Additionally, the \pkg{gamlss2} package introduces the capability to
  accommodate more than four parameter families. Essentially, this means that users can now
  specify any type of model using these new infrastructures, making the package highly
  flexible and accommodating to a wide range of modeling requirements.

  \packageIndices{gamlss2}
}

\author{
  \packageAuthor{gamlss2}

  Maintainer: \packageMaintainer{gamlss2}
}

\keyword{package}

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
  \code{\link{gamlss2}}, \code{\link{fake_formula}}
}

