\name{la}
\alias{la}
\alias{plot_lasso}

\title{Lasso Model Terms}

\description{
  Constructor function and plotting for Lasso penalized model terms for GAMLSS.
}

\usage{
## model term constructor function
la(x, type = 1, const = 1e-05, ...)

## plotting function
plot_lasso(x, terms = NULL,
  which = c("criterion", "coefficients"),
  zoom = c(3, 4), spar = TRUE, ...)
}

\arguments{
  \item{x}{For function \code{la()}, a numeric vector or matrix, or a formula.
    See the examples. For function \code{plot_lasso()}, an object returned from
    \code{\link{gamlss2}}.}
  \item{type}{Integer or character, the type of the Lasso penalty.
    \code{type = 1} or \code{type = "normal"} uses the normal penalty,
    \code{type = 2} or \code{type = "group"} the group penalty,
    \code{type = 3} or \code{type = "ordinal"} the ordinal fusion penalty and
    \code{type = 4} or \code{type = "nominal"} the nominal fusion penalty.}
  \item{const}{Numeric, the constant that is used for approximating the absolute function.}
  \item{terms}{Character or integer, the model term that should be plotted. The
    default \code{terms = NULL} is plotting all model terms.}
  \item{which}{Character, should the information criterion or the coefficient paths
    be plotted? See the examples.}
  \item{zoom}{Numeric vector of length 2, the zooming factors for plotting
    information criteria curves and coefficient paths. The first element sets
    the distance from the optimum shrinkage parameter lambda to the left side,
    and the second element to the right side, respectively.}
  \item{spar}{Logical, should plotting parameters be automatically set in
    \code{\link[graphics]{par}}?}
  \item{\dots}{For function \code{la()} further control arguments can be passed:
    The \code{criterion = "bic"} for shrinkage parameter selection, arguments
    for creating the \code{\link[stats]{model.matrix}} if the model term is specified
    using a formula, and \code{scale = TRUE}/\code{FALSE} to control the internal
    scaling of the design matrix. For function \code{plot_lasso()} arguments like
    \code{lwd}, \code{col}, \code{main}, etc., that control plotting parameters
    can be supplied. An additional ridge penalty (elastic net) can be added to
    each \code{la()} term be setting \code{add_ridge = TRUE}
    in the \code{\link{gamlss2}} call.}
}

\details{
  To implement the Lasso penalty, an approximation of the absolute value function is
  used, following the approach by Oelker and Tutz (2015). This enables the use of
  standard Newton-Raphson-type algorithms for estimation. Each Lasso model term has
  its own shrinkage parameter, allowing a mix of different penalty types within the
  model. The framework builds on the methodology of Groll et al. (2019),
  where coefficients are updated through iteratively reweighted least squares (IWLS).
  This is feasible due to the absolute function approximation, which results in
  a quadratic penalty matrix similar to that used in penalized splines. By default,
  the shrinkage parameters are selected using the Bayesian Information Criterion (BIC).

  The choice of \code{type} determines the structure that is penalized:
  \itemize{
    \item \code{"normal"} penalizes individual coefficients and is suitable for
      numeric covariates or general design matrices.
    \item \code{"group"} penalizes one complete block of coefficients jointly and
      is mainly useful for factor effects or other grouped design matrices.
    \item \code{"ordinal"} uses an adjacent-differences fusion penalty and is intended
      for ordered factor effects or matrices whose columns represent ordered levels.
    \item \code{"nominal"} uses a pairwise-differences fusion penalty and is intended
      for unordered factor effects or matrices whose columns represent unordered levels.
  }

  For numerical stability and comparable shrinkage across terms, \code{la()} can
  scale the underlying design matrix internally. The exact scaling depends on the
  selected penalty type. For \code{"ordinal"} and \code{"nominal"} terms, an
  additional ridge component is included by default.
}

\references{
  Andreas Groll, Julien Hambuckers, Thomas Kneib, and Nikolaus Umlauf (2019). Lasso-type penalization in
  the framework of generalized additive models for location, scale and shape.
  \emph{Computational Statistics & Data Analysis}.
  \doi{10.1016/j.csda.2019.06.005}

  Oelker Margreth-Ruth and Tutz Gerhard (2015). A uniform framework for combination of
  penalties in generalized structured models. \emph{Adv Data Anal Classif}.
  \doi{10.1007/s11634-015-0205-y}
}

\value{
  The \code{la()} function is used internally within \code{\link{gamlss2}}
  and provides the necessary details for estimating Lasso-type model terms.
  Essentially, it serves as a special model term, as outlined
  in \code{\link{specials}}.

  Currently, the \code{plot_lasso()} function does not return any output.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{specials}}.
}

\examples{\dontrun{## transform numeric to factor variables
air <- subset(airquality, !is.na(Ozone))

air$Tc <- cut(air$Temp, breaks = seq(50, 100, by = 10),
  include.lowest = TRUE)
air$Wc <- cut(air$Wind, breaks = seq(2, 22, by = 4),
  include.lowest = TRUE)

## set up the model formula for a BCT model
f <- Ozone ~ la(Tc,type=3) + la(Wc,type=3) | . | . | .

## estimation
m1 <- gamlss2(f, data = air, family = BCT)

## summary, shows the estimated degrees of freedom
## for each model term
summary(m1)

## plot estimated coefficients
plot(m1)

## plot information criteria curves
## for each model term
plot_lasso(m1)

## plot parameter paths
plot_lasso(m1, which = "coefficients")

## plot a single model term
plot_lasso(m1, which = "coefficients", term = 2)

## same with
plot_lasso(m1, which = "coefficients", term = "mu.la(Wc")

## zoom out
plot_lasso(m1, which = "coefficients", term = 2,
  zoom = c(8, 7))

## set names
plot_lasso(m1, which = "coefficients", term = 2,
  zoom = c(8, 7), names = c("A", "B", "C", "D"))

## set title
plot_lasso(m1, which = "coefficients", term = 2,
  zoom = c(8, 7), main = "Fused Lasso")

## simulated example using the nominal fused lasso
## and a matrix as argument for la()
set.seed(1328)

## number of observations and covariates
n <- 500
k <- 50

## model matrix
X <- matrix(rnorm(n * k), n, k)
colnames(X) <- paste0("x", 1:k)

## true coefficients
beta <- list(
  "mu" = rbinom(k, 1, 0.1),
  "sigma" = rbinom(k, 1, 0.1) * 0.3
)

## parameters
mu <- X \%*\% beta$mu
sigma <- exp(-1 + X \%*\% beta$sigma)

## response
y <- rnorm(n, mean = mu, sd = sigma)

## model formula with nominal fused lasso
f <- y ~ la(X,type=4) | la(X,type=4)

## estimate model incl. extra ridge penalty
## for all la() model terms
m2 <- gamlss2(f, add_ridge = TRUE)

## plot information criteria curves
plot_lasso(m2)

## coefficient paths
plot_lasso(m2, which = "coefficients")

## zoom out
plot_lasso(m2, which = "coefficients",
  zoom = c(8, 9))

## extract coefficients
cb <- coef(m2, full = TRUE)

## compare (without intercept)
cb_mu <- cb[grep("mu.", names(cb))][-1]
cb_sigma <- cb[grep("sigma.", names(cb))][-1]

## true positive rate
tp <- mean(c(cb_mu[beta$mu > 0] > 0,
  cb_sigma[beta$sigma > 0] > 0))

## false positive rate, needs threshold
thres <- 0.01
fp <- mean(c(abs(cb_mu[beta$mu == 0]) > thres,
  abs(cb_sigma[beta$sigma == 0]) > thres))
}
}

\keyword{regression}
\keyword{models}
\keyword{distribution}

