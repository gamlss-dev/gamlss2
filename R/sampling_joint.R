## Gaussian coefficient simulation from the cached joint precision.
sampling <- function(object, R = 100, antithetic = TRUE, ...)
{
  dots <- list(...)
  method <- dots$method %||% "joint"
  full <- isTRUE(dots$full)
  info <- joint_information(object, method = method)
  if(any(info$map$reason == "aliased", na.rm = TRUE))
    stop("Gaussian coefficient draws are unavailable with aliased coefficients")

  p <- info$dimension
  if(antithetic) {
    R0 <- ceiling(R / 2)
    Z <- matrix(rnorm(R0 * p), nrow = p, ncol = R0)
    Z <- cbind(Z, -Z)[, seq_len(R), drop = FALSE]
  } else {
    Z <- matrix(rnorm(R * p), nrow = p, ncol = R)
  }
  if(is.null(info$factor) || info$rank != p)
    stop("Gaussian coefficient draws require full-rank positive definite information")
  draws <- backsolve(info$factor, Z)
  draws <- sweep(draws, 1L, rowMeans(draws), "-")

  complete <- matrix(info$coefficients, nrow = length(info$coefficients),
    ncol = R)
  active <- which(info$map$active)
  complete[active, ] <- sweep(draws, 1L, info$active.coefficients, "+")
  keep <- if(full) seq_len(nrow(info$map)) else which(info$map$type == "linear")
  complete <- t(complete[keep, , drop = FALSE])
  colnames(complete) <- info$map$name[keep]
  complete
}
