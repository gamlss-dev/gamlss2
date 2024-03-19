## Function to compute multivariate normal samples
## given the maximum likelihood estimator.
sampling <- function(object, R = 100, ...)
{
  V <- vcov(object, ...)
  Cv <- chol(V)
  cb <- coef(object, ...)
  sc <- rnorm(R * length(cb))
  sc <- t(cb + t(Cv) %*% matrix(sc, nrow = length(cb), ncol = R))
  d <- drop(cb - apply(sc, 2, mean))
  sc <- t(t(sc) + d)
  colnames(sc) <- names(cb)
  return(sc)
}

