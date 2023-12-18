distfit <- function(..., plot = TRUE)
{
  b <- gamlss2(...)
  y <- b$y
  nd <- data.frame(seq(min(y), max(y), length = 100))
  yn <- as.character(formula(b$fake_formula, rhs = 0))[2L]
  names(nd) <- yn
  par <- predict(b, newdata = nd, type = "parameter")
  nd$dy <- family(b)$d(nd[[yn]], par)
  attr(y, "fitted.values") <- nd
  class(y) <- "distfit"
  if(plot) {
    plot(y)
  }
  return(y)
}

plot.distfit <- function(x, ...)
{
  h <- hist(x, plot = FALSE)
  d <- attr(x, "fitted.values")
  ylim <- range(h$density, d$dy)
  hist(x, breaks = "Scott", freq = FALSE, main = "",
    xlab = names(d)[1])
  lines(d, col = 4, lwd = 2)
}

