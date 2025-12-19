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
  y = NULL, what = "mu", breaks = seq(0, 1, by = 0.1),
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
  \item{newdata}{Optional data frame for out-of-sample calibration.
    If supplied, predictions and (optionally) the response \code{y} are
    evaluated on \code{newdata}.}
  \item{y}{Optional numeric or factor vector with the binary response
    (coded as 0/1 or a two-level factor). If omitted, the response is
    extracted from the first model.}
  \item{what}{Character, which parameter should be used for prediction.
    Typically \code{"mu"} for binary models.}
  \item{breaks}{Numeric vector of break points used to construct bins
    for the predicted probabilities. The range of \code{breaks} must
    cover the interval \code{[0, 1]}.}
  \item{minn}{Integer, the minimum number of observations required in
    each bin. Bins with fewer observations are dropped from the
    calibration curve.}
  \item{main}{Character, main title for the plot.}
  \item{xlab}{Character, label for the x-axis (predicted probabilities).}
  \item{ylab}{Character, label for the y-axis (observed proportions).}
  \item{plot}{Logical, should a calibration plot be produced?
    If \code{FALSE}, only the aggregated calibration data are returned.}
  \item{add_loess}{Logical, should a loess smooth be added for each model?}
  \item{smooth_n}{Integer, number of evaluation points for the loess curve.}
  \item{col}{Either \code{NULL}, a function, or a vector of colors.
    If \code{NULL}, distinct colors are generated via
    \code{colorspace::qualitative_hcl()}. If a function, it is called
    with the number of models to generate colors. Otherwise, a vector of
    colors is recycled to the number of models.}
  \item{lty}{Integer vector of line types used for the loess curves.
    If \code{NULL}, different line types are used for each model and
    recycled if necessary.}
  \item{legend}{Logical, should a legend be added when more than one
    model is supplied?}
  \item{pos}{Character, the position of the legend passed to
    \code{\link[graphics]{legend}}.}
  \item{xlim}{The x limits of the plot.}
  \item{ylim}{The y limits of the plot.}
}

\details{
  For each fitted model, predicted probabilities for \eqn{P(Y = 1 | X)}
  are obtained via \code{\link{predict}} with
  \code{type = "parameter"} and the requested \code{model} argument, if supplied.
  If the corresponding \code{\link{gamlss2.family}} object provides a
  \code{probabilities} method, this is used to transform the parameter
  vector into class probabilities. For families with more than two
  outcome categories, the first column is taken as the probability of
  the event.

  The predicted probabilities are then grouped into bins defined by
  \code{breaks}. Within each bin, the mean predicted probability and the
  observed proportion of 1s are computed, together with the number of
  observations in the bin. Bins with fewer than \code{minn}
  observations are dropped. The resulting points are plotted with size
  proportional to the square root of the bin size, and a reference line
  with slope 1 and intercept 0 is added for perfect calibration.
  Optionally, a loess smooth is fitted to the binned data for each model
  and added as a calibration curve.

  If multiple models are supplied, the binned calibration points and
  loess curves of all models are shown in the same plot, with different
  colors and line types. This allows a direct visual comparison of
  calibration across competing \code{\link{gamlss2}} specifications.
}

\value{
  Invisibly returns a data frame with the following columns:
  \itemize{
    \item{\code{interval}}{Factor with the probability intervals.}
    \item{\code{probs}}{Mean predicted probability in each bin.}
    \item{\code{y}}{Observed proportion of 1s in each bin.}
    \item{\code{n}}{Number of observations in each bin.}
    \item{\code{model}}{Character, model label. For a single model this
      column is dropped from the returned data frame.}
  }
  If \code{plot = FALSE}, the same data frame is returned without
  producing a plot.
}

\references{
  Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW (2019).
  Calibration: the Achilles heel of predictive analytics.
  \emph{BMC Medicine}, \bold{17}, 230.
  \doi{10.1186/s12916-019-1466-7}

  Steyerberg EW (2019).
  \emph{Clinical Prediction Models: A Practical Approach to Development,
  Validation, and Updating}. 2nd edition. Springer, New York.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{predict}}, \code{\link[graphics]{legend}}.
}

\examples{\dontrun{## simulate simple binary data
## and compare two models
set.seed(123)
n <- 1000
x1 <- runif(n, -3, 3)
x2 <- runif(n, -3, 3)
eta <- sin(x1) + 0.5 * x2
p <- plogis(eta)
y <- rbinom(n, size = 1, prob = p)
d <- data.frame(y = y, x1 = x1, x2 = x2)

## fit three competing models
m1 <- gamlss2(y ~ x1, family = BI, data = d)
m2 <- gamlss2(y ~ s(x1), family = BI, data = d)
m3 <- gamlss2(y ~ s(x1) + x2, family = BI, data = d)

## calibration on the training data
calibration(m1, m2, m3)

## extract calibration data without plotting
cal_tab <- calibration(m1, m2, m3, plot = FALSE)
head(cal_tab)
}}

