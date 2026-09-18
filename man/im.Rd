\name{im}
\alias{im}
\alias{special_fit.im}
\alias{special_predict.im.fitted}

\title{Convolutional Image Model Terms for GAMLSS}

\description{
  Constructs a convolutional neural-network image model term for use with
  \code{gamlss2}. The term learns convolutional filters, applies a ReLU
  activation and spatial-pyramid pooling, and fits a penalized linear output
  layer to the current working response.
}

\usage{
im(x, dim = NULL, ...)
}

\arguments{
  \item{x}{A single image variable. Images can be supplied as an \eqn{n} by \eqn{p}
    matrix of flattened images, an \eqn{n} by height by width (optionally by
    channel) array, or a list of equally sized image matrices or arrays.}
  \item{dim}{Optional integer vector \code{c(height, width)} for grayscale images or
    \code{c(height, width, channels)} for multi-channel images. It is required
    for flattened matrix input and can otherwise be inferred from the data.}
  \item{filters}{Integer number of convolutional filters. The default is 4.}
  \item{kernel}{One integer or two integers giving the convolutional kernel
    height and width. The default is \code{c(3, 3)}.}
  \item{pool}{Positive integer spatial-pyramid pooling levels. The default is
    \code{c(1, 2, 4)}, truncated to the image dimensions.}
  \item{epochs}{Number of optimization epochs for the first backfitting update.
    The default is 8.}
  \item{warm_epochs}{Number of optimization epochs for subsequent updates.
    The default is 3.}
  \item{learning_rate}{Positive Adam learning rate. The default is 0.003.}
  \item{decay}{Non-negative ridge penalty for the convolutional and output
    parameters. The default is 0.001.}
  \item{beta1,beta2}{Adam exponential-decay parameters in \code{[0, 1)}.
    Their defaults are 0.9 and 0.999.}
  \item{eps}{Positive numerical stabilization constant for Adam. The default
    is \code{1e-8}.}
  \item{clip}{Optional positive bound for gradient clipping. The default is 5.}
  \item{scale}{Logical; should each image channel be centered and scaled using
    the training images? The default is \code{TRUE}.}
  \item{seed}{Optional finite numeric seed used when initializing the filters.}
  \item{trace}{Logical; should the CNN optimization progress be printed? The
    default is \code{FALSE}.}
  \item{...}{Additional named CNN controls listed above. Unknown controls are
    rejected.}
}

\details{
  The image term is intended for images stored as one observation per row. For
  a flattened matrix, use \code{dim} to specify the image dimensions; for an
  array or image list, the dimensions are checked and retained. The same image
  representation and dimensions must be supplied when predicting on new data.

  The model applies a same-padded convolution, a rectified-linear activation,
  and spatial-pyramid pooling. The pooled features are mapped to the current
  working response by a penalized linear output layer. Convolutional parameters
  are initialized on the first update and warm-started on later backfitting
  updates. The fitted contribution is centered before it is added to the
  distributional predictor.

  This term is a CNN-based alternative to \code{\link{si}}, which represents
  images through a penalized tensor-product coefficient surface. The fitted
  model can be inspected with \code{\link{plot.gamlss2}} and used for new
  images with \code{\link{predict.gamlss2}}.
}

\value{
  An object of classes \code{"special"} and \code{"im"}. During fitting, the
  corresponding \code{special_fit.im()} method returns the fitted CNN,
  centered fitted values, transfer information for warm starts, and an
  effective degrees-of-freedom estimate.
}

\seealso{
  \code{\link{si}}, \code{\link{special_terms}}, \code{\link{gamlss2}}
}

\examples{\dontrun{
set.seed(1328)
n <- 120
nr <- 16
nc <- 16

## images with a three-by-three pattern at random locations
images <- array(rnorm(n * nr * nc, sd = 0.25), c(n, nr, nc))
locations <- data.frame(
  row = sample(2:(nr - 2), n, replace = TRUE),
  col = sample(2:(nc - 2), n, replace = TRUE)
)
amplitude <- runif(n, 0.5, 1.5)
for(i in seq_len(n)) {
  images[i, locations$row[i]:(locations$row[i] + 2),
    locations$col[i]:(locations$col[i] + 2)] <-
    images[i, locations$row[i]:(locations$row[i] + 2),
      locations$col[i]:(locations$col[i] + 2)] + amplitude[i]
}

## known spatial response surface
truth <- outer(seq_len(nr), seq_len(nc), function(row, col)
  exp(-((row - 5) / 3)^2 - ((col - 11) / 4)^2))
truth <- truth / max(truth)
d <- data.frame(y = vapply(seq_len(n), function(i) {
  amplitude[i] * mean(truth[locations$row[i]:(locations$row[i] + 2),
    locations$col[i]:(locations$col[i] + 2)])
}, numeric(1)) + rnorm(n, sd = 0.2))
image(images[1, , ], col = hcl.colors(20, "Inferno"),
  main = "Example image")
d$image <- I(matrix(images, nrow = n))

m <- gamlss2(
  y ~ im(image, dim = c(nr, nc), filters = 3, epochs = 15,
    warm_epochs = 2, seed = 1904) | 1,
  data = d,
  family = NO
)

plot(m)
predict(m, newdata = d[1:5, , drop = FALSE])

## predicted response map for a moving image pattern
map.grid <- expand.grid(row = 2:(nr - 2), col = 2:(nc - 2))
map.images <- array(0, c(nrow(map.grid), nr, nc))
for(i in seq_len(nrow(map.grid))) {
  map.images[i, map.grid$row[i]:(map.grid$row[i] + 2),
    map.grid$col[i]:(map.grid$col[i] + 2)] <- 1
}
map.data <- data.frame(y = rep(0, nrow(map.grid)))
map.data$image <- I(matrix(map.images, nrow = nrow(map.grid)))
map.pred <- predict(m, parameter = "mu", newdata = map.data)
true.map <- vapply(seq_len(nrow(map.grid)), function(i) {
  sum(truth[map.grid$row[i]:(map.grid$row[i] + 2),
    map.grid$col[i]:(map.grid$col[i] + 2)]) / 9
}, numeric(1))
true.map <- matrix(true.map,
  nrow = length(unique(map.grid$row)))
pred.map <- matrix(map.pred,
  nrow = length(unique(map.grid$row)))

op <- par(mfrow = c(1, 2))
image(unique(map.grid$row), unique(map.grid$col), true.map,
  col = hcl.colors(20, "Inferno"),
  xlab = "image row", ylab = "image column",
  main = "True response map")
image(unique(map.grid$row), unique(map.grid$col), pred.map,
  col = hcl.colors(20, "Inferno"),
  zlim = range(true.map, pred.map),
  xlab = "image row", ylab = "image column",
  main = "Predicted response map")
par(op)
}}

\keyword{models}
\keyword{neural}
