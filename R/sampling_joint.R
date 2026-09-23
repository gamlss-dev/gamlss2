## Gaussian coefficient simulation from the full covariance.
sampling <- function(object, R = 100, antithetic = TRUE, ...)
{
  dots <- list(...)
  method <- dots$method %||% "joint"
  method <- match.arg(method, c("joint", "working", "numeric"))
  full <- isTRUE(dots$full)
  info <- vcov_information(object, method = method)
  draws <- vcov_draws(info, R, antithetic = antithetic, center = TRUE)
  blocks <- info$blocks
  linear <- unlist(lapply(blocks, function(z) {
    k <- which(vapply(z, function(x) identical(x$type, "linear"), logical(1L)))
    if(!length(k)) integer() else z[[k[1L]]]$index
  }), use.names = FALSE)
  keep <- if(full) seq_len(info$dimension) else linear
  answer <- t(draws[keep, , drop = FALSE])
  colnames(answer) <- names(info$coefficients)[keep]
  answer
}
