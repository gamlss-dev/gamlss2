## Fast joint covariance conditional on the fitted smoothing parameters.
vcov.gamlss2 <- function(object,
  type = c("vcov", "cor", "se", "coef"), full = FALSE,
  method = c("joint", "working", "numeric"), ...)
{
  type <- match.arg(type)
  method <- match.arg(method)
  if(type == "coef") {
    return(coef(object, full = full, dropall = FALSE))
  }

  info <- joint_information(object, method = method, ...)
  keep <- if(full) seq_len(nrow(info$map)) else which(info$map$type == "linear")
  Vfull <- expand_vcov(info)
  V <- Vfull[keep, keep, drop = FALSE]
  if(type == "cor") {
    se <- sqrt(diag(V))
    V <- V / outer(se, se)
    diag(V)[is.finite(se) & se > 0] <- 1
  }
  if(type == "se") {
    d <- diag(V)
    scale <- if(any(is.finite(d))) max(1, abs(d[is.finite(d)])) else 1
    tolerance <- sqrt(.Machine$double.eps) * scale
    negative <- !is.na(d) & d < -tolerance
    if(any(negative)) {
      warning("materially negative coefficient variance; standard error is NA",
        call. = FALSE)
      d[negative] <- NA_real_
    }
    d[!is.na(d) & d < 0] <- 0
    V <- sqrt(d)
    names(V) <- rownames(Vfull)[keep]
  }
  V
}
