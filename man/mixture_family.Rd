\name{mixture_family}
\alias{mixture_family}

\title{Construct a Finite Mixture Distribution Family}

\description{
  Constructs a \code{"gamlss2.family"} representing a finite mixture of
  two or more component distribution families. Component distribution
  parameters and the mixing probabilities may all be modeled using
  separate additive predictors.
}

\usage{
mixture_family(families = NO, k = NULL, reference = 1L,
  prefix = "c", initialize = TRUE)
}

\arguments{
  \item{families}{
    A single distribution family or a list of component families. Each
    entry may be a family-generating function, a family object, or a
    character family name accepted by the \pkg{gamlss2} family
    infrastructure. If a single family is supplied, it is repeated
    \code{k} times.
  }
  \item{k}{
    Number of mixture components. For a single family, \code{NULL}
    defaults to two components. For a list of families, \code{k} is
    inferred from its length and, if supplied, must agree with that
    length.
  }
  \item{reference}{
    Integer identifying the reference component for the mixing-weight
    parameterization. Its unnormalized odds are fixed at one.
  }
  \item{prefix}{
    Either one non-empty character prefix, from which component labels
    such as \code{"c1"} and \code{"c2"} are formed, or one distinct
    label per component. Component parameter names are formed by
    joining the label and original parameter name with a dot.
  }
  \item{initialize}{
    Logical. If \code{TRUE}, location parameters named \code{mu},
    \code{mean}, or \code{location} receive separated empirical-quantile
    starting values. Other starting-value functions are inherited from
    their component families.
  }
}

\details{
  For component densities \eqn{f_j(y)} and probabilities \eqn{\pi_j},
  the mixture density is
  \deqn{f(y) = \sum_{j=1}^k \pi_j f_j(y).}
  It is evaluated using a log-sum-exp calculation to avoid numerical
  underflow.

  Component parameters are ordered component by component. For example,
  a two-component normal mixture has parameters
  \code{c1.mu}, \code{c1.sigma}, \code{c2.mu}, \code{c2.sigma}, and
  \code{pi2}. Consequently, a \code{\link{gamlss2}} formula for this
  family has five right-hand sides in that order.

  Mixing probabilities use a baseline-logit parameterization analogous
  to \code{\link{MN}}. For every non-reference component \eqn{j}, the
  parameter \code{pij} has a log link and represents positive odds
  \eqn{w_j} relative to the reference component. With the reference
  odds fixed at one,
  \deqn{\pi_j = \frac{w_j}{1 + \sum_{\ell \ne r} w_\ell}, \qquad \pi_r = \frac{1}{1 + \sum_{\ell \ne r} w_\ell}.}
  The \code{probabilities()} family function converts these odds to
  normalized probabilities.

  Analytical first derivatives are assembled from the component-family
  scores. If
  \deqn{r_{ij} = \frac{\pi_{ij}f_j(y_i)}{\sum_\ell \pi_{i\ell}f_\ell(y_i)}}
  denotes a posterior responsibility, the score for a component
  parameter is its component score multiplied by \eqn{r_{ij}}. The
  mixing-predictor score is \eqn{r_{ij} - \pi_{ij}}. Missing Hessians
  are supplied numerically by the standard family completion
  infrastructure.

  The returned family includes \code{cdf()}, \code{quantile()},
  \code{random()}, \code{mean()}, and \code{variance()} whenever all
  component families provide the corresponding operation. All components
  must use the same distribution type, for example all continuous or all
  discrete. Mixture quantiles are obtained by numerically inverting the
  mixture CDF.

  As with finite mixtures generally, the likelihood can have local
  optima and equivalent solutions obtained by exchanging component
  labels. Separated initialization helps but does not impose an ordering
  constraint. Multiple fits with different starting values can still be
  useful for difficult data.
}

\value{
  An object of class \code{"gamlss2.family"}. In addition to the standard
  family elements, it contains:
  \itemize{
    \item \code{probabilities(par, n = NULL)}: normalized prior mixing
      probabilities.
    \item \code{responsibilities(par, y)}: posterior component
      probabilities conditional on observations.
    \item \code{components}: the completed component families.
    \item \code{component_map}: mappings from original component
      parameter names to mixture parameter names.
  }
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2.family}}
}

\examples{
\dontrun{## waiting times in the faithful data are distinctly bimodal
data("faithful", package = "datasets")

## sort for plotting
d <- faithful[order(faithful$eruptions), ]

## construct and estimate a two-normal mixture
fam <- mixture_family(NO, k = 2)
m <- gamlss2(waiting ~ s(eruptions) | . | . | . | ., data = d, family = fam)

## estimated effects
plot(m)

## predict parameters
par <- predict(m)
head(par)

## prior and posterior component probabilities
prior <- family(m)$probabilities(par)
head(prior)
posterior <- family(m)$responsibilities(par, d$waiting)
head(posterior)

## classify observations by their largest
## posterior probability and show estimated quantiles.
fit <- quantile(m)

pc <- max.col(as.matrix(posterior))
ccols <- c(2, 4)

plot(
  waiting ~ eruptions,
  data = d,
  col = ccols[pc],
  pch = 19,
  xlab = "Eruption duration (minutes)",
  ylab = "Waiting time (minutes)",
  ylim = range(d$waiting, fit)
)

matlines(d$eruptions, fit, col = 1, lty = 1)

legend(
  "topleft", paste("Component", 1:2),
  col = ccols, pch = 19, bty = "n"
)

## mixture regression with a numeric and a binary predictor.
data("mtcars", package = "datasets")
cars <- transform(
  mtcars,
  am = factor(am, levels = 0:1, labels = c("automatic", "manual"))
)

## both component means depend on weight and transmission,
## component scales and the mixing probability are constant
m <- gamlss2(
  mpg ~ s(wt) + am | 1 | s(wt) + am | 1 | 1,
  data = cars,
  family = fam
)

cars_par <- predict(m)
cars_posterior <- family(m)$responsibilities(cars_par, cars$mpg)
cars_component <- max.col(as.matrix(cars_posterior))
table(cars$am, cars_component)

## plot estimated classes and component mean regressions,
## point shape and line type distinguish automatic
## from manual transmission.
wt_grid <- seq(min(cars$wt), max(cars$wt), length.out = 100)
newcars <- expand.grid(wt = wt_grid, am = levels(cars$am))
newcars$am <- factor(newcars$am, levels = levels(cars$am))
cars_fit <- predict(m, newdata = newcars)

plot(
  mpg ~ wt,
  data = cars,
  col = ccols[cars_component],
  pch = c(1, 16)[cars$am],
  xlab = "Weight (1000 lbs)",
  ylab = "Miles per gallon",
  ylim = range(cars$mpg, cars_fit)
)

for(j in 1:2) {
  for(a in seq_along(levels(cars$am))) {
    keep <- newcars$am == levels(cars$am)[a]
    lines(
      newcars$wt[keep], cars_fit[keep, paste0("c", j, ".mu")],
      col = ccols[j], lty = a, lwd = 2
    )
  }
}

legend(
  "topright", paste("Component", 1:2),
  col = ccols, lty = 1, lwd = 2, bty = "n"
)
legend(
  "bottomleft", levels(cars$am),
  col = 1, pch = c(1, 16), lty = 1:2, bty = "n"
)
}}

\keyword{distribution}
\keyword{mixture}
\keyword{regression}

