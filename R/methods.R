## Mean function.
mean.gamlss2 <- function(x, ...) {
  stopifnot(requireNamespace("distributions3"))
  d <- distributions3::prodist(x, ...)
  mean(d)
}

## Median function.
median.gamlss2 <- function(x, ...) {
  quantile(x, ..., probs = 0.5)
}

## Variance function.
variance.gamlss2 <- function(x, ...) {
  stopifnot(requireNamespace("distributions3"))
  d <- distributions3::prodist(x, ...)
  variance(d)
}

## Skewness function.
skewness.gamlss2 <- function(x, ...) {
  stopifnot(requireNamespace("distributions3"))
  d <- distributions3::prodist(x, ...)
  skewness(d)
}

## Kurtosis function.
kurtosis.gamlss2 <- function(x, ...) {
  stopifnot(requireNamespace("distributions3"))
  d <- distributions3::prodist(x, ...)
  kurtosis(d)
}

## Density function.
pdf.gamlss2 <- function(d, x, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  if(missing(x))
    x <- model.response(model.frame(d))
  d <- distributions3::prodist(d, ...)
  pdf(d, x, drop = drop, ...)
}

## Log-density function.
log_pdf.gamlss2 <- function(d, x, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  if(missing(x))
    x <- model.response(model.frame(d))
  d <- distributions3::prodist(d, ...)
  pdf(d, x, drop = drop, log = TRUE)
}

## Cumulative distribution function.
cdf.gamlss2 <- function(d, x, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  if(missing(x))
    x <- model.response(model.frame(d))
  d <- distributions3::prodist(d, ...)
  cdf(d, x, drop = drop, ...)
}

## Random numbers.
random.gamlss2 <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  d <- distributions3::prodist(x, ...)
  random(d, n, drop = drop, ...)
}

is_discrete.gamlss2 <- function(d, ...) {
  stopifnot(requireNamespace("distributions3"))
  f <- family(d)
  if(is.null(f$type)) stop(sprintf("the type is not implemented for the %s family",  f$family[1L]))
  return(tolower(f$type) == "discrete")
}

is_continuous.gamlss2 <- function(d, ...) {
  stopifnot(requireNamespace("distributions3"))
  f <- family(d)
  if(is.null(f$type)) stop(sprintf("the type is not implemented for the %s family", f$family[1L]))
  return(tolower(f$type) == "continuous")
}

support.gamlss2 <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  d <- distributions3::prodist(d, ...)
  s <- quantile(d, probs = c(0, 1), elementwise = FALSE)
  distributions3::make_support(s[, 1L], s[, 2L], d, drop = drop)
}

## Quantiles.
quantile.gamlss2 <- function(x,
  probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
  variable = NULL, newdata = NULL, plot = FALSE, data = TRUE, n = 100L, ...)
{
  mf <- model.frame(x)
  rn <- all.vars(formula(x$fake_formula, rhs = 0))
  mf <- mf[, !(names(mf) %in% rn), drop = FALSE]
  if(is.null(variable)) {
    par <- predict(x, newdata = if(is.null(newdata)) mf else newdata, type = "parameter", drop = FALSE)
  } else {
    if(length(av <- all.vars(formula(x$fake_formula, lhs = 0))) > 1L)
      stop("variable plot can only be generated for single covariate models!")
    if(length(variable) > 1L)
      stop("argument variable is specified wrong!")
    if(is.logical(variable))
      variable <- av
    if(!is.character(variable))
      variable <- names(mf)[variable]
    variable <- grep(variable, names(mf), value = TRUE, fixed = TRUE)[1L]
    xr <- range(mf[, variable])
    if(is.null(newdata)) {
      newdata <- list()
      for(j in names(mf))
        newdata[[j]] <- rep(mean(mf[[j]]), lenghth = n)
      newdata[[variable]] <- seq(xr[1L], xr[2L], length = n)
      newdata <- as.data.frame(newdata)
    } else {
      newdata <- newdata[order(newdata[[variable]]), ]
    }
    par <- predict(x, newdata = newdata, type = "parameter", drop = FALSE)
  }
  f <- NULL
  for(p in probs)
    f <- cbind(f, family(x)$quantile(p, par))
  colnames(f) <- paste0(probs * 100, "%")
  if(plot) {
    if(is.null(variable)) {
      ind <- 1:nrow(f)
      ind <- ind / max(ind) * 100
      plot(range(c(x$y, f)), range(ind), type = "n",
        xlab = rn, ylab = "%")
      if(data) {
        points(sort(x$y), 1:length(x$y) / length(x$y) * 100,
          col = rgb(0.1, 0.1, 0.1, alpha = 0.3), pch = 16)
      }
      col <- list(...)$col
      if(is.null(col))
        col <- rev(colorspace::heat_hcl(ncol(f) + 1L))[-1L]
      col <- rep(col, length.out = ncol(f))
      for(j in 1:ncol(f)) {
        lines(ind ~ sort(f[, j]), col = col[j], lwd = 2)
      }
      legend <- list(...)$legend
      if(is.null(legend))
        legend <- TRUE
      if(isTRUE(legend)) {
        legend("topleft", colnames(f), col = col, lwd = 2, bty = "n")
      }
    } else {
      if(data) {
        plot(mf[[variable]], x$y,
          col = rgb(0.1, 0.1, 0.1, alpha = 0.3), pch = 16,
          xlab = variable, ylab = rn)
      }
      col <- list(...)$col
      if(is.null(col))
        col <- rev(colorspace::heat_hcl(ncol(f) + 1L))[-1L]
      col <- rep(col, length.out = ncol(f))
      matplot(newdata[[variable]], f, type = "l", col = col, lwd = 2, lty = 1,
        add = data, xlab = if(data) "" else variable,
        ylab = if(data) "" else rn)
      legend <- list(...)$legend
      if(is.null(legend))
        legend <- TRUE
      if(isTRUE(legend)) {
        legend("topleft", colnames(f), col = col, lwd = 2, bty = "n")
      }
    }
  }
  cn <- colnames(f)
  if(!is.null(variable)) {
    f <- as.data.frame(cbind(newdata[[variable]], f))
    colnames(f) <- c(x, cn)
  } else {
    f <- as.data.frame(f)
    colnames(f) <- cn
    if(ncol(f) < 2L)
      f <- as.numeric(f[, 1])
  }
  if(plot) {
    return(invisible(f))
  } else {
    return(f)
  }
}

