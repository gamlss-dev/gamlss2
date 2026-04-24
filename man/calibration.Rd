\name{calibration}
\alias{calibration}

\title{Calibration Plots for Binary Responses}

\description{
  Compute and plot calibration curves for models with 0/1 responses.
  The function can handle one or several fitted \code{\link{gamlss2}} models
  and compares their calibration in a single plot.
}

\usage{
calibration(..., newdata = NULL,
  y = NULL, parameter = NULL, breaks = seq(0, 1, by = 0.1),
  minn = 20, main = "Calibration plot",
  xlab = "Predicted probability",
  ylab = "Observed proportion", plot = TRUE,
  add_loess = TRUE, smooth_n = 200,
  col = NULL, lty = NULL, legend = TRUE, pos = "topleft",
  xlim = NULL, ylim = NULL)
}

\arguments{
  \item{\dots}{One or several fitted \code{\link{gamlss2}} model objects
    to be assessed.}
  \item{newdata}{Optional data frame for out-of-sample calibration.}
  \item{y}{Optional numeric or factor vector with the binary response.}
  \item{parameter}{Character, which parameter should be used for prediction.}
  \item{breaks}{Numeric vector of break points used to construct bins
    for the predicted probabilities.}
  \item{minn}{Integer, the minimum number of observations required in
    each bin.}
  \item{main}{Character, main title for the plot.}
  \item{xlab}{Character, label for the x-axis.}
  \item{ylab}{Character, label for the y-axis.}
  \item{plot}{Logical, should a calibration plot be produced?}
  \item{add_loess}{Logical, should a loess smooth be added for each model?}
  \item{smooth_n}{Integer, number of evaluation points for the loess curve.}
  \item{col}{Colors used for the models.}
  \item{lty}{Line types used for the loess curves.}
  \item{legend}{Logical, should a legend be added when more than one
    model is supplied?}
  \item{pos}{Character, legend position passed to
    \code{\link[graphics]{legend}}.}
  \item{xlim}{The x limits of the plot.}
  \item{ylim}{The y limits of the plot.}
}

\details{
  For each fitted model, predicted probabilities are obtained via
  \code{\link{predict}} using \code{type = "parameter"} and the selected
  \code{parameter}. If the corresponding \code{\link{gamlss2.family}} object
  provides a \code{probabilities()} method, this is used to transform the
  parameter vector into probabilities. For multi-column outputs, the first
  column is used.

  Predicted probabilities are grouped into bins defined by \code{breaks}.
  Within each bin, the mean predicted probability, the observed proportion
  of 1s, and the number of observations are computed. Bins with fewer than
  \code{minn} observations are dropped.
}

\value{
  Invisibly returns a data frame with the columns \code{interval},
  \code{probs}, \code{y}, \code{n}, and \code{model}. For a single model,
  the \code{model} column is dropped.
}

\references{
  Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW (2019).
  Calibration: the Achilles heel of predictive analytics.
  \emph{BMC Medicine}, \bold{17}, 230.
  \doi{10.1186/s12916-019-1466-7}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{predict}}
}

\examples{\dontrun{set.seed(123)
n <- 400
x1 <- runif(n, -3, 3)
x2 <- runif(n, -3, 3)
eta <- sin(x1) + 0.5 * x2
p <- plogis(eta)
y <- rbinom(n, size = 1, prob = p)
d <- data.frame(y = y, x1 = x1, x2 = x2)

m1 <- gamlss2(y ~ x1, family = BI, data = d)
m2 <- gamlss2(y ~ s(x1), family = BI, data = d)
m3 <- gamlss2(y ~ s(x1) + x2, family = BI, data = d)

calibration(m1, m2, m3)
head(calibration(m1, m2, m3, plot = FALSE))
}}
