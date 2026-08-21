\name{gamlss2-package}
\alias{gamlss2-package}

\docType{package}

\title{GAMLSS Modeling with Advanced Flexible Infrastructures}

\description{
  Next generation infrastructure for generalized additive models for
  location, scale, and shape (GAMLSS) and distributional regression
  more generally. The package provides a fresh reimplementation of
  the classic 'gamlss' package while being more modular and
  facilitating the creation of advanced terms and models.
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
  Nikolaus Umlauf,
  Mikis Stasinopoulos,
  Robert Rigby,
  Achim Zeileis,
  Reto Stauffer,
  Gillian Heller,
  Fernanda de Bastiani,
  Andreas Mayr.

  Maintainer: Nikolaus Umlauf \email{Nikolaus.Umlauf@uibk.ac.at}.
}

\keyword{package}

\references{
To cite gamlss2 in publications use:

  Umlauf N, Stasinopoulos M, De Bastiani F, Heller G, Kneib T, Mayr A,
  Rigby R, Stauffer R, Zeileis A (2025). \dQuote{gamlss2: The Next Generation
  of Distributional Regression.} In O'Sullivan D, Sweeney J, Moghaddam
  S, Burke K (eds.), \emph{Proceedings of the 39th International Workshop on
  Statistical Modelling}, 312--317. \url{https://iwsm2025.ie/}

For the underlying methods cite:

  Rigby RA, Stasinopoulos DM (2005). \dQuote{Generalized Additive Models for Location,
  Scale and Shape (with Discussion).} \emph{Journal of the Royal Statistical Society,
  Series C (Applied Statistics)}, \bold{54}, 507--554.
  \doi{10.1111/j.1467-9876.2005.00510.x}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{fake_formula}}
}
