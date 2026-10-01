## Gaussian coefficient simulation from the full covariance.
sampling <- function(object, R = 100, antithetic = TRUE, ...,
  unconditional = FALSE, sandwich = FALSE)
{
  dots <- list(...)
  method <- dots$method %||% "joint"
  method <- match.arg(method, c("joint", "working", "numeric"))
  full <- isTRUE(dots$full)
  info <- vcov_information(object, method = method,
    unconditional = unconditional, sandwich = sandwich,
    .details = sandwich)
  if(sandwich && length(info$details$rho)) {
    se <- sqrt(diag(info$details$smoothing))
    if(any(!is.finite(se)) ||
        any(2 * qnorm(0.975) * se > log(100)))
      stop(paste0("RS sandwich Gaussian draws are unreliable: a 95% ",
        "log-smoothing-parameter interval spans more than a factor of ",
        "100. The local covariance does not determine predictive tail ",
        "probabilities for this fit; full RS refits are needed to ",
        "estimate such bands."), call. = FALSE)
  }
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
