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
  \pkg{gamlss2} fits an additive predictor for each parameter of a response
  distribution, using linear terms, smooth terms from \pkg{mgcv}, and
  user-defined model terms (see \code{\link{specials}}). Response distributions
  are specified through \code{\link{gamlss2.family}} objects,
  including families with more than four parameters. Fitting algorithms can
  be selected through \code{\link{gamlss2_control}}. Fitted models provide
  methods for summaries, prediction, plotting, and diagnostics. Predictions
  can be returned as \pkg{distributions3} objects for evaluating fitted
  distributions (see \code{\link{predict.gamlss2}}).

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
  of Distributional Regression}. In O'Sullivan D, Sweeney J, Moghaddam
  S, Burke K (eds.), \emph{Proceedings of the 39th International Workshop on
  Statistical Modelling}, 312--317. Department of Mathematics and Statistics,
  University of Limerick, Limerick, Ireland.
  \url{https://iwsm2025.ie/wp-content/uploads/2025/07/IWSM2025_Limerick_Proceedings.pdf}

For the underlying methods cite:

  Rigby RA, Stasinopoulos DM (2005). \dQuote{Generalized Additive Models for Location,
  Scale and Shape (with Discussion)}. \emph{Journal of the Royal Statistical Society,
  Series C (Applied Statistics)}, \bold{54}(3), 507--554.
  \doi{10.1111/j.1467-9876.2005.00510.x}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{fake_formula}}
}
