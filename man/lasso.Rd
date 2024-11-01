\name{la}
\alias{la}
\alias{plot_lasso}

\title{Lasso Model Terms}

\description{
  Constructor function and plotting for Lasso penalized model terms for GAMLSS.
}

\usage{
## Model term constructor function.
la(x, type = 1, const = 1e-05, ...)

## Plotting function.
plot_lasso(x, terms = NULL,
  which = c("criterion", "coefficients"),
  zoom = c(3, 4), spar = TRUE, ...)
}

\arguments{
  \item{x}{For function \code{la()}, a numeric vector or matrix, or a formula.
    See the examples. For function \code{plot_lasso()}, an object returned from
    \code{\link{gamlss2}}.
  \item{type}{Integer or character, the type of the Lasso penalty.
    \code{type = 1} or \code{type = "normal") uses the normal penalty,
    \code{type = 2} or \code{type = "group") the group penalty,
    \code{type = 3} or \code{type = "nominal") the nominal fusion penalty and
    \code{type = 4} or \code{type = "ordinal") the ordinal fusion penalty.
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
    using a formula. For function \code{plot_lasso()} arguments like
    \code{lwd}, \code{col}, \code{main}, etc., that control plotting parameters
    can be supplied.}
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
  \code{\link{gamlss2}}, \code{link{specials}}.
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
\dontrun{## load data
data("rent", package = "gamlss.data")

## transform numeric to factor variables
rent$Flc <- cut(rent$Fl, breaks = seq(20, 160, by = 10),
  include.lowest = TRUE)
rent$Ac <- cut(rent$A, breaks = seq(1890, 1990, by = 10),
  include.lowest = TRUE)

## set up the model formula for a BCT model
f <- R ~ la(Flc,type=4) + la(Ac,type=4) + la(loc,type=4) |
  la(Flc,type=4) + la(Ac,type=4) + la(loc,type=4) |
  la(Flc,type=4) + la(Ac,type=4) + la(loc,type=4)

## estimation
b <- gamlss2(f, data = rent, family = BCT)

## summary, shows the estimated degrees of freedom
## for each model term
summary(b)

## plot estimated coefficients
plot(b)

## plot information criteria curves
## for each model term.
plot_lasso(b)

## plot parameter paths.
plot_lasso(b, which = "coefficients")

## plot a single model term.
plot_lasso(b, which = "coefficients", term = 5)

## same with
plot_lasso(b, which = "coefficients", term = "sigma.la(Ac")

## zoom out
plot_lasso(b, which = "coefficients", term = 5,
  zoom = c(8, 7))

## set names
plot_lasso(b, which = "coefficients", term = 5,
  zoom = c(8, 7), names = c("A", "B", "C"))

## set title
plot_lasso(b, which = "coefficients", term = 5,
  zoom = c(8, 7), main = "Fused Lasso")

## simulated example using the normal lasso
## and a matrix as argument for la()
set.seed(123)

## number of observations and covariates
n <- 500
k <- 50

## model matrix and coefficients
X <- matrix(rnorm(n * k), n, k)
colnames(X) <- paste0("x", 1:k)
beta <- cbind("mu" = rbinom(k, 1, 0.1), "sigma" = rbinom(k, 1, 0.1) * 0.3)

## parameters
mu <- X %*% beta[, "mu"]
sigma <- exp(1 + X %*% beta[, "sigma"])

## response
y <- rnorm(n, mean = mu, sd = sigma)

## model formula
f <- y ~ la(X) | la(X)

## estimate model
b <- gamlss2(f)

## plot information criteria curves
plot_lasso(b)

## coefficient paths
plot_lasso(b, which = "coefficients")

## zoom out
plot_lasso(b, which = "coefficients",
  zoom = c(8, 9))

## extract coefficients
cb <- coef(b, full = TRUE)

## compare (without intercept)
cb_mu <- cb[grep("mu.", names(cb))][-1]
cb_sigma <- cb[grep("sigma.", names(cb))][-1]

print(cbind(beta[, "mu"], round(cb_mu, 2)))
print(cbind(beta[, "sigma"], round(cb_sigma, 2)))
}
}

\keyword{regression}
\keyword{models}
\keyword{distribution}

