\name{MN}
\alias{MN}

\title{Multinomial Logit Family}

\description{
The \code{MN()} family implements a multinomial logit model for
categorical responses with \eqn{k \ge 2} unordered categories.
It can be used within \code{\link{gamlss2}} to model category
probabilities via separate additive predictors for each
non-reference category (first level per default).
}

\usage{
MN(k)
}

\arguments{
  \item{k}{An integer specifying the number of response categories.}
}

\details{
  The response must be a factor with exactly \code{k} levels.
  The first level is treated as the reference category.

  Let \eqn{Y \in \{1, \ldots, k\}} denote the response and let
  \eqn{\eta_j} be the additive predictor for category \eqn{j}
  (\eqn{j = 2, \ldots, k}). The model is defined as

  \eqn{\Pr(Y = 1) = \frac{1}{1 + \sum_{j=2}^k \exp(\eta_j)}},

  \eqn{\Pr(Y = j) = \frac{\exp(\eta_j)}{1 + \sum_{m=2}^k \exp(\eta_m)}, \qquad j = 2, \ldots, k}.

  Thus, each predictor models the log-odds relative to the
  reference category:

  \eqn{\log \frac{\Pr(Y = j)}{\Pr(Y = 1)} = \eta_j}.

  The link function for all non-reference categories is the log link,
  i.e., the internal parameterization uses positive weights
  \eqn{\exp(\eta_j)}.

  The family provides analytical first and second derivatives of the
  log-likelihood with respect to the predictors. The Hessian entries
  correspond to the negative second derivatives, which are required by
  the Rigby-Stasinopoulos (\code{\link{RS}}) backfitting algorithm. In particular,

  \itemize{
    \item Diagonal elements:
      \eqn{-\partial^2 \ell / \partial \eta_j^2 = p_j (1 - p_j)}
    \item Off-diagonal elements:
      \eqn{-\partial^2 \ell / \partial \eta_j \partial \eta_m = - p_j p_m}
  }

  where \eqn{p_j = \Pr(Y = j)}.

  The family object provides a \code{probabilities()} method returning
  the estimated category probabilities for all observations, see the example below.
}

\value{
  A \code{"gamlss2.family"} object to be used with \code{\link{gamlss2}}.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2.family}}
}

\examples{\dontrun{data("marital.nz", package = "VGAM")

## estimate model, same predictor
## for each category using | .
m <- gamlss2(mstatus ~ s(age) | . | ., data = marital.nz, family = MN(4))

## plot estimated effects
plot(m)

## predict probabilities
par <- predict(m)
probs <- family(m)$probabilities(par)
names(probs) <- levels(marital.nz$mstatus)

## plot fitted probabilities
i <- order(marital.nz$age)
matplot(marital.nz$age[i], probs[i, ], type = "l", lty = 1, lwd = 2,
  xlab = "Age", ylab = "Probabilities")
legend("center", names(probs), lty = 1, lwd = 2, col = 1:4)

## plot quantile residuals
plot(m, which = c("hist-resid", "qq-resid", "wp-resid"))
}}

\keyword{regression}
\keyword{categorical}
\keyword{multinomial}
\keyword{distribution}

