\name{si}
\alias{si}
\alias{smooth.construct.si.smooth.spec}
\alias{Predict.matrix.si.effect}

\title{Smooth Image Model Terms for GAMLSS}

\description{
  Constructs a smooth image model term for \code{gamlss2}. The image
  coefficient surface is represented by a tensor product of marginal
  B-spline bases with difference penalties.
}

\usage{
si(x, dim = NULL, k = c(8, 8), degree = c(3, 3),
  m = c(2, 2), sp = NULL, ...)
}

\arguments{
  \item{x}{
    A matrix-valued data-frame column containing one flattened image per row.
    The number of columns must equal the product of \code{dim}.
  }
  \item{dim}{
    Integer vector \code{c(height, width)} for grayscale images or
    \code{c(height, width, channels)} for multi-channel images.
  }
  \item{k}{
    One or two integers giving the marginal B-spline basis dimensions.
    The default is \code{c(8, 8)}.
  }
  \item{degree}{
    One or two non-negative integers giving the marginal B-spline degrees.
    The default is cubic splines, \code{c(3, 3)}.
  }
  \item{m}{
    One or two difference-penalty orders, each smaller than the corresponding
    value of \code{k}. The default is \code{c(2, 2)}.
  }
  \item{sp}{
    Optional vector of one or two non-negative smoothing parameters. If
    supplied, the two marginal smoothing parameters are kept fixed.
  }
  \item{\dots}{
    Further controls passed to the smooth specification.
  }
}

\details{
  The image contribution is a linear functional of a tensor-product
  coefficient surface. For image \eqn{x_i(h,w)} and coefficient surface
  \eqn{\beta(h,w)}, the fitted contribution is
  \deqn{f_i = \sum_h \sum_w x_i(h,w) \beta(h,w).}
  The surface is expanded in marginal B-spline bases, and the image
  values are projected onto that basis without explicitly constructing
  the full image-by-coefficient Kronecker matrix.

  For multi-channel images, one coefficient surface is fitted per channel
  with common marginal smoothing parameters. The fitted term is centered
  before it is added to the distributional predictor.

  The compact data representation is usually created with
  \code{dat$image <- I(matrix(images, nrow = n))}, where \code{images} is
  an \code{n} by height by width (and optionally channels) array.
}

\value{
  The constructor returns a smooth specification object of class
  \code{"si.smooth.spec"}. After construction by \pkg{mgcv},
  the corresponding effect contains the tensor-product design and
  penalty matrices used for fitting and prediction.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{special_terms}},
  \code{\link[mgcv]{smooth.construct}}
}

\examples{\dontrun{
set.seed(1328)
n <- 365*3
nr <- 30
nc <- 40
lat <- seq(40, 55, length.out = nr)
lon <- seq(-10, 20, length.out = nc)

## artificial daily temperature fields
clim <- outer(lat, lon,
  function(lat, lon) 20 - 0.35 * (lat - 40) + 0.04 * lon)
pw <- outer(lat, lon,
  function(lat, lon) exp(-((lat - 48) / 4)^2 - ((lon + 3) / 6)^2))
pe <- outer(lat, lon,
  function(lat, lon) exp(-((lat - 46) / 4)^2 - ((lon - 14) / 5)^2))
pn <- outer(as.numeric(scale(lat)), rep(1, nc))

temp <- array(NA_real_, c(n, nr, nc))
for(i in seq_len(n)) {
  temp[i, , ] <- clim +
    rnorm(1, sd = 2.5) +
    rnorm(1, sd = 3.0) * pw +
    rnorm(1, sd = 3.0) * pe +
    rnorm(1, sd = 1.0) * pn +
    matrix(rnorm(nr * nc, sd = 0.25), nr, nc)
}

## ozone also depends on humidity, wind,
## and two temperature regions
hum <- runif(n, 35, 95)
wind <- rgamma(n, shape = 4, scale = 1.2)
tw <- apply(temp[, 4:7, 2:5, drop = FALSE], 1L, mean)
te <- apply(temp[, 3:6, 10:13, drop = FALSE], 1L, mean)
f_im <- 1.8 * (te - mean(te)) - 1.2 * (tw - mean(tw))
f_hum <- -0.006 * (hum - 55)^2
f_wind <- 7 * exp(-0.5 * ((wind - 3) / 2)^2)

d <- data.frame(
  ozone = 45 + f_hum + f_wind + f_im + rnorm(n, sd = 2),
  hum = hum,
  wind = wind
)
d$temp <- I(matrix(temp, nrow = n))

## three temperature fields
zlim <- range(temp[1:3, , ])
op <- par(mfrow = c(1, 3))
for(i in 1:3) {
  image(lon, lat, t(temp[i, , ]),
    zlim = zlim, asp = 1,
    col = hcl.colors(30, "Inferno", rev = TRUE),
    xlab = "lon", ylab = "lat",
    main = paste("Temperature, day", i))
}
par(op)

m <- gamlss2(
  ozone ~ s(hum) + s(wind) +
    si(temp, dim = c(nr, nc), k = 10),
  data = d,
  family = NO
)

## plot the fitted temperature-to-ozone image effect
ss <- m$specials[["si(temp)"]]
cf <- m$fitted.specials$mu[["si(temp)"]]$coefficients
## each pixel is evaluated as a one-unit image perturbation,
## this visualizes the fitted coefficient surface,
## not a reconstructed field
P <- kronecker(ss$Bw, ss$Bh)
q <- attr(ss, "qrc")
P <- t(qr.qty(q, t(P))[(attr(ss, "nCons") + 1L):ncol(P), , drop = FALSE])
eff <- drop(P \%*\% cf)
image(lon, lat, t(matrix(eff, nr, nc)),
  col = hcl.colors(30, "Inferno", rev = TRUE),
  xlab = "lon", ylab = "lat", main = "Fitted si() image effect")

## predicted distributions and quantiles
nd <- d[1:5, , drop = FALSE]
predict(m, newdata = nd)
quantile(m, newdata = nd)
}}

\keyword{models}
\keyword{smooth}
