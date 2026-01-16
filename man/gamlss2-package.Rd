\name{gamlss2-package}
\alias{gamlss2-package}

\docType{package}

\title{GAMLSS Modeling with Advanced Flexible Infrastructures}

\description{
  Next generation infrastructure for generalized additive models for
  location, scale, and shape (GAMLSS) and general distributional regression
  models. The package provides a fresh reimplementaton of
  the classic 'gamlss' package while being more modular and
  facilitating the creation of advanced terms and models.
}

\details{
  The primary purpose of this package is to facilitate the creation of advanced infrastructure
  designed to enhance the Generalized Additive Models for Location Scale and Shape (GAMLSS,
  Rigby and Stasinopoulos 2005) modeling framework. Notably, the \pkg{gamlss2} package represents a
  significant overhaul of its predecessor, \pkg{gamlss}, with a key emphasis on improving estimation
  speed and incorporating more adaptable infrastructures. These enhancements enable the seamless
  integration of various estimation algorithms into GAMLSS, including gradient boosting, Bayesian estimation,
  regression trees, and forests, fostering a more versatile and powerful modeling environment.

  Moreover, the package expands its compatibility by supporting all model terms from the base
  R \pkg{mgcv} package. Additionally, the \pkg{gamlss2} package introduces the capability to
  accommodate more than four parameter families. Essentially, this means that users can now
  specify any type of model using these new infrastructures, making the package highly
  flexible and accommodating to a wide range of modeling requirements.
}

\author{
  Nikolaus Umlauf [aut, cre] \email{Nikolaus.Umlauf@uibk.ac.at}
  (\href{https://orcid.org/0000-0003-2160-9803}{ORCID: 0000-0003-2160-9803}) \cr
  Mikis Stasinopoulos [fnd, cph] \email{d.stasinopoulos@gre.ac.uk}
  (\href{https://orcid.org/0000-0003-2407-5704}{ORCID: 0000-0003-2407-5704}) \cr
  Robert Rigby [fnd, cph] \email{r.rigby@gre.ac.uk}
  (\href{https://orcid.org/0000-0003-3853-1707}{ORCID: 0000-0003-3853-1707}) \cr
  Achim Zeileis [ctb] \email{Achim.Zeileis@R-project.org}
  (\href{https://orcid.org/0000-0003-0918-3766}{ORCID: 0000-0003-0918-3766}) \cr
  Reto Stauffer [ctb] \email{Reto.Stauffer@uibk.ac.at}
  (\href{https://orcid.org/0000-0002-3798-5507}{ORCID: 0000-0002-3798-5507}) \cr
  Gillian Heller [ctb] \email{Gillian.Heller@sydney.edu.au}
  (\href{https://orcid.org/0000-0003-1270-1499}{ORCID: 0000-0003-1270-1499}) \cr
  Fernanda de Bastiani [ctb] \email{fernandadebastiani@gmail.com}
  (\href{https://orcid.org/0000-0001-8532-639X}{ORCID: 0000-0001-8532-639X}) \cr
  Andreas Mayr [ctb] \email{Andreas.Mayr@uni-marburg.de}
  (\href{https://orcid.org/0000-0001-7106-9732}{ORCID: 0000-0001-7106-9732})
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
  \code{\link{gamlss2}}
}

