\name{im}
\alias{im}
\alias{special_fit.im}
\alias{special_predict.im.fitted}
\title{Interaction Model Terms}
\description{Constructs an interaction model term for use with \code{gamlss2}.}
\usage{im(x, dim = NULL, ...)}
\arguments{\item{x}{A model term specifying the variables entering the interaction.}\item{dim}{Optional dimension of the interaction representation.}\item{...}{Additional arguments passed to the interaction-term implementation.}}
\details{The term is fitted and predicted through the special-term backfitting interface. See \code{\link{special_terms}} for the common fitting and prediction contract.}
\value{An object representing an interaction special term.}
\seealso{
  \code{\link{special_terms}}
}

\keyword{models}
