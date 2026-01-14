\name{elm}
\alias{elm}

\title{Extreme Learning Machine Model Terms}

\description{
  Constructor function for Extreme Learning Machine (ELM) model terms for GAMLSS.
}

\usage{
## Model term constructor function.
elm(x, k = 50, a = "tanh", ...)
}

\arguments{
  \item{x}{A numeric vector or matrix, a factor, or a formula.
    If \code{x} is a formula, a design matrix is created using
    \code{\link[stats]{model.matrix}}. See the examples.}
  \item{k}{Integer, number of hidden units (random features) used to build the
    ELM design matrix.}
  \item{a}{Character, activation function used for the hidden units.
    Supported options are \code{"logistic"}, \code{"tanh"} (default),
    \code{"relu"} and \code{"identity"}.}
  \item{\dots}{Further control arguments can be passed:
    \code{criterion = "bic"} (default) for shrinkage parameter selection and
    \code{scale = TRUE} (default) for internal scaling of the design matrix.
    Further arguments are passed to \code{\link[stats]{model.matrix}} if
    \code{x} is specified using a formula.}
}

\details{
  The ELM model term constructs a randomized single-hidden-layer representation
  of the covariate(s) in \code{x}. Internally, a design matrix \code{Z} is built
  from \code{x} (including an intercept column). Random weights are sampled and
  transformed through an activation function to obtain the hidden-layer design
  matrix \code{X}. The columns of \code{X} are centered before estimation.

  For numerical stability and comparability across terms, the design matrix
  \code{Z} can be scaled internally (\code{scale = TRUE}). If \code{x} is a factor,
  a QR-based group scaling is applied; otherwise, a column-wise normal scaling is
  used. The applied scaling is stored internally and automatically reused during
  prediction.

  The returned object is a special model term used within \code{\link{gamlss2}},
  following the interface described in \code{\link{specials}}.
}

\value{
  The \code{elm()} function is used internally within \code{\link{gamlss2}} and
  provides the necessary details for estimating an ELM model term. Essentially,
  it serves as a special model term, as outlined in \code{\link{specials}}.
}

\references{
  Huang GB, Zhu QY, Siew CK (2006).
  Extreme Learning Machine: Theory and Applications.
  \emph{Neurocomputing}, \bold{70}(1--3), 489--501.
  \doi{10.1016/j.neucom.2005.12.126}

  Huang GB, Zhu QY, Siew CK (2004).
  Extreme Learning Machine: A New Learning Scheme of Feedforward Neural Networks.
  In: \emph{Proceedings of IJCNN 2004}, \bold{2}, 985--990.
  \doi{10.1109/IJCNN.2004.1380068}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{specials}}.
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
\dontrun{
## load data
data("SpirometryUS")

## subset for female
d <- subset(SpirometryUS, gender == "Female")

## note, inner weight are sampled, set
## the seed for reprodicibility
set.seed(1328)

## formula for all 4 parameters
f <- fev1 ~ elm(age) | . | . | .

## estimate model
m <- gamlss2(f, data = d, family = BCT)

## estimated effects
plot(m)

## predict quantiles
qu <- c(0.025, seq(0.1, 0.9, by = 0.1), 0.975)
fit <- quantile(m, probs = qu)

## plot
plot(fev1 ~ age, data = d)
i <- order(d$age)
matplot(d$age[i], fit[i, ],
  type = "l", lty = 1, lwd = 2,
  col = c(2, rep(4, ncol(fit) - 1)),
  add = TRUE)

## main effects and interactions
f <- fev1 ~ s(age) + s(height) + s(weight) +
  elm(~ age + height + weight, k = 200) | . | . | .

m <- gamlss2(f, data = d, family = BCT)

## summary to inspect effect of interactions
## of the elm()s
summary(m)

## plot main effects
plot(m)

## prediction is handled automatically via
## the special term interface
n <- 50
nd <- with(d, expand.grid(
  "age" = seq(min(age), max(age), length = n),
  "weight" = seq(min(weight), max(height), length = n)
))
nd$height <- mean(d$height)

## compute lower 2.5% limit
nd$fit <- quantile(m, newdata = nd, probs = 0.025)

## visualize
n <- length(unique(nd$age))
age <- sort(unique(nd$age))
weight <- sort(unique(nd$weight))

z <- matrix(nd$fit, n, n)

image(age, weight, z,
      col = hcl.colors(100, "YlOrRd"),
      xlab = "age", ylab = "weight")
contour(age, weight, z, add = TRUE)
}
}

