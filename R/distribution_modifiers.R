## Utilities for composing gamlss2 families without changing the base family.

.modifier_rows <- function(par, n = NULL)
{
  par <- as.list(par)
  lengths <- lengths(par)
  lengths <- lengths[lengths > 1L]
  nr <- if(length(lengths)) max(lengths) else 1L
  if(!is.null(n)) nr <- max(nr, n)
  nr
}

.modifier_recycle <- function(x, n)
{
  if(length(x) == n) x else rep(x, length.out = n)
}

.modifier_subset <- function(par, i, n)
{
  par <- as.list(par)
  lapply(par, function(x) {
    x <- .modifier_recycle(x, n)
    x[i]
  })
}

.modifier_call_cdf <- function(family, par, y, lower.tail = TRUE,
  log.p = FALSE, ...)
{
  cdf <- family$cdf
  if(!is.function(cdf))
    stop("the base family needs a $cdf() function", call. = FALSE)

  fml <- names(formals(cdf))
  args <- list(par = par, y = y)
  direct.tail <- "lower.tail" %in% fml || "..." %in% fml
  direct.log <- "log.p" %in% fml || "log" %in% fml

  if(direct.tail)
    args$lower.tail <- lower.tail
  if("log.p" %in% fml)
    args$log.p <- log.p && (lower.tail || direct.tail)
  else if("log" %in% fml)
    args$log <- log.p && (lower.tail || direct.tail)

  dots <- list(...)
  dots[c("par", "y", "lower.tail", "log.p", "log")] <- NULL
  value <- as.numeric(do.call(cdf, c(args, dots)))

  if(!direct.tail && !lower.tail) {
    if(direct.log && log.p) {
      args[[if("log.p" %in% fml) "log.p" else "log"]] <- FALSE
      value <- as.numeric(do.call(cdf, c(args, dots)))
    }
    value <- pmin(pmax(value, 0), 1)
    value <- if(log.p) log1p(-value) else 1 - value
  }

  value
}

.modifier_call_quantile <- function(family, par, p, lower.tail = TRUE,
  log.p = FALSE, ...)
{
  quantile <- family$quantile
  if(!is.function(quantile))
    stop("the base family needs a $quantile() function", call. = FALSE)

  fml <- names(formals(quantile))
  direct.tail <- "lower.tail" %in% fml || "..." %in% fml
  args <- list(par = par, p = p)
  if(direct.tail)
    args$lower.tail <- lower.tail
  if("log.p" %in% fml)
    args$log.p <- log.p
  else if("log" %in% fml)
    args$log <- log.p

  if(!direct.tail || (!"log.p" %in% fml && !"log" %in% fml)) {
    probability <- if(log.p) exp(p) else p
    if(!lower.tail) probability <- 1 - probability
    args$p <- probability
    args$lower.tail <- NULL
    args$log.p <- args$log <- NULL
  }
  dots <- list(...)
  dots[c("par", "p", "lower.tail", "log.p", "log")] <- NULL
  as.numeric(do.call(quantile, c(args, dots)))
}

.modifier_mass <- function(family, par, y, log = FALSE, ...)
{
  if(is.function(family$mass))
    return(family$mass(par = par, y = y, log = log, ...))

  if(!identical(tolower(family$type[1L]), "discrete")) {
    value <- rep.int(if(log) -Inf else 0, length(y))
    value[is.na(y)] <- NA_real_
    return(value)
  }

  family$pdf(par = par, y = y, log = log, ...)
}

.modifier_cdf_left <- function(family, par, y, log.p = FALSE, ...)
{
  if(is.function(family$cdf_left))
    return(family$cdf_left(par = par, y = y, log.p = log.p, ...))

  if(!identical(tolower(family$type[1L]), "discrete"))
    return(.modifier_call_cdf(family, par, y, log.p = log.p, ...))

  distribution <- .modifier_call_cdf(family, par, y, ...)
  mass <- .modifier_mass(family, par, y, ...)
  value <- pmax(distribution - mass, 0)
  if(log.p) log(value) else value
}

.modifier_logdiff <- function(a, b)
{
  value <- rep.int(-Inf, max(length(a), length(b)))
  a <- .modifier_recycle(a, length(value))
  b <- .modifier_recycle(b, length(value))
  ok <- is.finite(a) & (b == -Inf | b <= a)
  value[ok] <- a[ok] + log1p(-exp(b[ok] - a[ok]))
  value[a == b] <- -Inf
  value[a == 0 & b == -Inf] <- 0
  value
}

.modifier_log_interval <- function(family, par, lower, upper,
  lower.left = FALSE, upper.left = FALSE, ...)
{
  lower.log <- if(lower.left) {
    .modifier_cdf_left(family, par, lower, log.p = TRUE, ...)
  } else {
    .modifier_call_cdf(family, par, lower, log.p = TRUE, ...)
  }
  upper.log <- if(upper.left) {
    .modifier_cdf_left(family, par, upper, log.p = TRUE, ...)
  } else {
    .modifier_call_cdf(family, par, upper, log.p = TRUE, ...)
  }

  lower.log <- pmin(lower.log, 0)
  upper.log <- pmin(upper.log, 0)
  from.cdf <- .modifier_logdiff(upper.log, lower.log)
  lower.survival <- if(!lower.left ||
      !identical(tolower(family$type[1L]), "discrete")) {
    .modifier_call_cdf(
      family, par, lower, lower.tail = FALSE, log.p = TRUE, ...
    )
  } else {
    log1p(-exp(lower.log))
  }
  upper.survival <- if(!upper.left ||
      !identical(tolower(family$type[1L]), "discrete")) {
    .modifier_call_cdf(
      family, par, upper, lower.tail = FALSE, log.p = TRUE, ...
    )
  } else {
    log1p(-exp(upper.log))
  }
  from.survival <- .modifier_logdiff(lower.survival, upper.survival)

  use.survival <- upper.log > log(0.5) & is.finite(from.survival)
  from.cdf[use.survival] <- from.survival[use.survival]
  from.cdf
}

.modifier_links <- function(family)
{
  links <- lapply(family$links, make.link2)
  names(links) <- family$names
  links
}

.modifier_fd_score <- function(fun, par, parameter, link,
  step = .Machine$double.eps^(1 / 3))
{
  eta <- link$linkfun(par[[parameter]])
  plus <- minus <- par
  plus[[parameter]] <- link$linkinv(eta + step)
  minus[[parameter]] <- link$linkinv(eta - step)
  (fun(plus) - fun(minus)) / (2 * step)
}

.modifier_score_hessians <- function(score, family,
  step = .Machine$double.eps^(1 / 3))
{
  links <- .modifier_links(family)
  hessian <- list()
  for(first in family$names) {
    for(second in family$names) {
      if(match(first, family$names) > match(second, family$names))
        next
      key <- if(first == second) first else paste(first, second, sep = ":")
      reverse <- paste(second, first, sep = ":")
      hessian[[key]] <- local({
        perturb <- first
        evaluate <- score[[second]]
        link <- links[[first]]
        function(par, y, ...) {
          eta <- link$linkfun(par[[perturb]])
          plus <- minus <- par
          plus[[perturb]] <- link$linkinv(eta + step)
          minus[[perturb]] <- link$linkinv(eta - step)
          -(evaluate(par = plus, y = y, ...) -
            evaluate(par = minus, y = y, ...)) / (2 * step)
        }
      })
      if(first != second)
        hessian[[reverse]] <- hessian[[key]]
    }
  }
  hessian
}

.modifier_add_derivatives <- function(family, adjustment)
{
  links <- .modifier_links(family)
  score <- list()

  for(parameter in family$names) {
    base.score <- family$score[[parameter]]
    score[[parameter]] <- local({
      id <- parameter
      bscore <- base.score
      link <- links[[parameter]]
      function(par, y, ...) {
        fun <- function(p) adjustment(p, y)
        bscore(par = par, y = y, ...) +
          .modifier_fd_score(fun, par, id, link)
      }
    })
  }

  ## Some base families provide Fisher rather than observed information.
  ## Differentiating the modified score keeps the correction coherent.
  list(score = score, hessian = .modifier_score_hessians(score, family))
}

.modifier_initialize <- function(family, response)
{
  if(is.null(family$initialize))
    return(NULL)
  initialize <- lapply(family$initialize, function(fun) {
    if(!is.function(fun)) return(NULL)
    local({
      initialize <- fun
      function(y, ...) initialize(response(y), ...)
    })
  })
  names(initialize) <- names(family$initialize)
  initialize
}

.modifier_family <- function(base, name, type = base$type)
{
  list(
    family = name,
    names = base$names,
    links = base$links,
    type = type,
    map2par = base$map2par,
    base_family = base
  )
}

.modifier_no_moment <- function(label)
{
  force(label)
  function(...)
    stop(label, " is not available for this modified family", call. = FALSE)
}

## Construct the push-forward by a fixed monotone transformation.
transform_family <- function(family = NO, transform, inverse,
  log_jacobian, increasing = TRUE, support = c(-Inf, Inf),
  valid.response = NULL, name = NULL)
{
  family <- complete_family(family)
  if(!identical(tolower(family$type[1L]), "continuous"))
    stop("transform_family() currently requires a continuous base family",
      call. = FALSE)
  if(!is.function(transform) || !is.function(inverse) ||
      !is.function(log_jacobian)) {
    stop("'transform', 'inverse', and 'log_jacobian' must be functions",
      call. = FALSE)
  }
  if(!is.logical(increasing) || length(increasing) != 1L ||
      is.na(increasing))
    stop("'increasing' must be TRUE or FALSE", call. = FALSE)
  if(!is.numeric(support) || length(support) != 2L || anyNA(support) ||
      support[1L] >= support[2L])
    stop("'support' must contain increasing response bounds",
      call. = FALSE)
  if(is.null(valid.response))
    valid.response <- function(y) is.na(y) | (y >= support[1L] &
      y <= support[2L] & !is.na(suppressWarnings(inverse(y))))
  if(!is.function(valid.response))
    stop("'valid.response' must be a function", call. = FALSE)

  label <- if(is.null(name)) paste0("Transform(", family$family[1L], ")")
    else as.character(name)[1L]
  fam <- .modifier_family(family, label, "continuous")

  fam$pdf <- function(par, y, log = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy0 <- .modifier_recycle(y, n)
    x <- suppressWarnings(inverse(yy0))
    valid <- valid.response(yy0)
    value <- rep.int(-Inf, n)
    value[is.na(yy0)] <- NA_real_
    if(any(valid & !is.na(valid))) {
      i <- which(valid & !is.na(valid))
      pp <- .modifier_subset(par, i, n)
      yy <- x[i]
      jacobian <- .modifier_recycle(log_jacobian(
        yy0[i]
      ), length(i))
      value[i] <- family$pdf(par = pp, y = yy, log = TRUE, ...) + jacobian
    }
    if(log) value else exp(value)
  }

  fam$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    value <- numeric(n)
    inside <- !is.na(yy) & yy > support[1L] & yy < support[2L]
    if(any(inside)) {
      pp <- .modifier_subset(par, inside, n)
      x <- suppressWarnings(inverse(yy[inside]))
      value[inside] <- .modifier_call_cdf(
        family, pp, x, lower.tail = increasing, log.p = FALSE, ...
      )
    }
    value[!is.na(yy) & yy >= support[2L]] <- 1
    value[is.na(yy)] <- NA_real_
    if(!lower.tail) value <- 1 - value
    if(log.p) log(value) else value
  }

  fam$quantile <- function(par, p, lower.tail = TRUE, log.p = FALSE, ...) {
    probability <- if(log.p) exp(p) else p
    if(!lower.tail) probability <- 1 - probability
    if(!increasing) probability <- 1 - probability
    transform(family$quantile(par = par, p = probability, ...))
  }

  fam$random <- function(par, n, ...)
    transform(family$random(par = par, n = n, ...))
  fam$valid.response <- function(y) all(is.na(y) | valid.response(y))
  fam$initialize <- .modifier_initialize(family, inverse)
  fam$score <- lapply(family$score, function(fun) {
    local({
      score <- fun
      function(par, y, ...) score(par = par, y = inverse(y), ...)
    })
  })
  fam$hessian <- lapply(family$hessian, function(fun) {
    local({
      hessian <- fun
      function(par, y, ...) hessian(par = par, y = inverse(y), ...)
    })
  })
  fam$update <- if(is.function(family$update)) {
    function(par, y, eta, which)
      family$update(par = par, y = inverse(y), eta = eta, which = which)
  } else NULL
  fam$log_likelihood <- function(par, y, ...)
    sum(fam$pdf(par = par, y = y, log = TRUE, ...), na.rm = TRUE)
  fam$rqres <- NULL
  fam$support <- NULL
  fam$mean <- fam$variance <- NULL
  class(fam) <- "gamlss2.family"
  fam
}

## Logistic push-forward of a continuous family.
logit_family <- function(family = NO)
{
  transform_family(
    family = family,
    transform = stats::plogis,
    inverse = stats::qlogis,
    log_jacobian = function(y) -log(y) - log1p(-y),
    support = c(0, 1),
    valid.response = function(y) is.na(y) | (is.finite(y) & y > 0 & y < 1),
    name = paste0("Logit(", complete_family(family)$family[1L], ")")
  )
}

## Construct a family conditioned to a fixed interval.
truncate_family <- function(family = NO, lower = -Inf, upper = Inf,
  include = c(FALSE, TRUE))
{
  family <- complete_family(family)
  if(length(lower) != 1L || length(upper) != 1L || is.na(lower) ||
      is.na(upper) || lower >= upper)
    stop("'lower' and 'upper' must define a nonempty fixed interval",
      call. = FALSE)
  if(!is.logical(include) || length(include) != 2L || anyNA(include))
    stop("'include' must contain two logical values", call. = FALSE)

  lower.left <- isTRUE(include[1L])
  upper.left <- !isTRUE(include[2L])
  inside <- function(y) {
    lo <- if(include[1L]) y >= lower else y > lower
    hi <- if(include[2L]) y <= upper else y < upper
    lo & hi
  }
  log.normalizer <- function(par, n) {
    .modifier_log_interval(
      family, par, rep.int(lower, n), rep.int(upper, n),
      lower.left = lower.left, upper.left = upper.left
    )
  }
  adjustment <- function(par, y)
    -log.normalizer(par, max(length(y), .modifier_rows(par)))

  fam <- .modifier_family(
    family,
    paste0("Truncated(", family$family[1L], ")")
  )
  fam$pdf <- function(par, y, log = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    value <- family$pdf(par = par, y = yy, log = TRUE, ...) -
      log.normalizer(par, n)
    value[!inside(yy)] <- -Inf
    if(log) value else exp(value)
  }
  fam$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    numerator <- .modifier_log_interval(
      family, par, rep.int(lower, n), yy,
      lower.left = lower.left, upper.left = FALSE, ...
    )
    value <- exp(numerator - log.normalizer(par, n))
    value[yy < lower | (yy == lower & !include[1L])] <- 0
    value[yy >= upper] <- 1
    value <- pmin(pmax(value, 0), 1)
    if(!lower.tail) value <- 1 - value
    if(log.p) log(value) else value
  }
  fam$quantile <- function(par, p, lower.tail = TRUE, log.p = FALSE, ...) {
    probability <- if(log.p) exp(p) else p
    if(!lower.tail) probability <- 1 - probability
    n <- max(length(probability), .modifier_rows(par))
    probability <- .modifier_recycle(probability, n)
    a <- if(lower.left) {
      .modifier_cdf_left(family, par, rep.int(lower, n))
    } else {
      .modifier_call_cdf(family, par, rep.int(lower, n))
    }
    log.z <- log.normalizer(par, n)
    value <- rep.int(NA_real_, n)
    use.survival <- a > 0.5
    if(any(!use.survival)) {
      i <- !use.survival
      target <- a[i] + probability[i] * exp(log.z[i])
      value[i] <- .modifier_call_quantile(
        family, .modifier_subset(par, i, n),
        pmin(pmax(target, 0), 1), ...
      )
    }
    if(any(use.survival)) {
      i <- use.survival
      pp <- .modifier_subset(par, i, n)
      log.survival <- if(lower.left &&
          identical(tolower(family$type[1L]), "discrete")) {
        left <- .modifier_cdf_left(
          family, pp, rep.int(lower, sum(i)), log.p = TRUE
        )
        log1p(-exp(left))
      } else {
        .modifier_call_cdf(
          family, pp, rep.int(lower, sum(i)),
          lower.tail = FALSE, log.p = TRUE
        )
      }
      removed <- log(probability[i]) + log.z[i]
      target <- .modifier_logdiff(log.survival, removed)
      value[i] <- .modifier_call_quantile(
        family, pp, exp(target), lower.tail = FALSE, ...
      )
    }

    ## Some legacy quantile wrappers clip tail probabilities.  Fall back to
    ## vectorized inversion only for candidates that fail their CDF check.
    if(identical(tolower(family$type[1L]), "continuous")) {
      fitted.probability <- fam$cdf(par, value)
      bad <- !is.na(probability) & probability > 0 & probability < 1 &
        (!is.finite(value) | abs(fitted.probability - probability) > 1e-7)
      if(any(bad)) {
        pp <- .modifier_subset(par, bad, n)
        target <- probability[bad]
        candidate <- value[bad]
        candidate[!is.finite(candidate)] <- 0
        left <- if(is.finite(lower)) rep.int(lower, sum(bad)) else
          candidate - 1
        right <- if(is.finite(upper)) rep.int(upper, sum(bad)) else
          candidate + 1
        width <- rep.int(1, sum(bad))
        if(!is.finite(lower)) {
          for(iteration in seq_len(60L)) {
            active <- fam$cdf(pp, left) >= target
            if(!any(active)) break
            left[active] <- left[active] - width[active]
            width[active] <- 2 * width[active]
          }
        }
        width[] <- 1
        if(!is.finite(upper)) {
          for(iteration in seq_len(60L)) {
            active <- fam$cdf(pp, right) < target
            if(!any(active)) break
            right[active] <- right[active] + width[active]
            width[active] <- 2 * width[active]
          }
        }
        for(iteration in seq_len(80L)) {
          midpoint <- (left + right) / 2
          move <- fam$cdf(pp, midpoint) < target
          left[move] <- midpoint[move]
          right[!move] <- midpoint[!move]
        }
        value[bad] <- (left + right) / 2
      }
    }
    value
  }
  fam$random <- function(par, n, ...)
    fam$quantile(par = par, p = stats::runif(n), ...)
  fam$mass <- function(par, y, log = FALSE, ...) {
    value <- .modifier_mass(family, par, y, log = TRUE, ...) -
      log.normalizer(par, max(length(y), .modifier_rows(par)))
    value[!inside(y)] <- -Inf
    if(log) value else exp(value)
  }
  fam$cdf_left <- function(par, y, log.p = FALSE, ...) {
    value <- fam$cdf(par, y, ...)
    mass <- fam$mass(par, y, ...)
    value <- pmax(value - mass, 0)
    if(log.p) log(value) else value
  }
  derivatives <- .modifier_add_derivatives(family, adjustment)
  fam$score <- derivatives$score
  fam$hessian <- derivatives$hessian
  fam$update <- NULL
  fam$log_likelihood <- function(par, y, ...)
    sum(fam$pdf(par = par, y = y, log = TRUE, ...), na.rm = TRUE)
  fam$valid.response <- function(y) all(is.na(y) | inside(y))
  fam$rqres <- NULL
  fam$support <- NULL
  fam$mean <- fam$variance <- NULL
  class(fam) <- "gamlss2.family"
  fam
}

## Construct the observed law after fixed-boundary or right censoring.
censor_family <- function(family = NO, lower = -Inf, upper = Inf)
{
  family <- complete_family(family)
  if(length(lower) != 1L || length(upper) != 1L || is.na(lower) ||
      is.na(upper) || lower >= upper)
    stop("'lower' and 'upper' must define a nonempty fixed interval",
      call. = FALSE)
  if(!is.finite(lower) && !is.finite(upper))
    return(.right_censor_family(family))

  lower.logprob <- function(par, n) {
    if(!is.finite(lower)) return(rep.int(-Inf, n))
    .modifier_call_cdf(
      family, par, rep.int(lower, n), log.p = TRUE
    )
  }
  upper.logprob <- function(par, n) {
    if(!is.finite(upper)) return(rep.int(-Inf, n))
    left <- .modifier_cdf_left(
      family, par, rep.int(upper, n), log.p = TRUE
    )
    log1p(-exp(pmin(left, 0)))
  }

  fam <- .modifier_family(
    family,
    paste0("Censored(", family$family[1L], ")"),
    if(identical(family$type, "continuous")) "mixed" else family$type
  )
  fam$pdf <- function(par, y, log = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    value <- family$pdf(par = par, y = yy, log = TRUE, ...)
    if(is.finite(lower)) value[yy == lower] <- lower.logprob(par, n)[yy == lower]
    if(is.finite(upper)) value[yy == upper] <- upper.logprob(par, n)[yy == upper]
    value[yy < lower | yy > upper] <- -Inf
    if(log) value else exp(value)
  }
  fam$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
    yy <- .modifier_recycle(y, max(length(y), .modifier_rows(par)))
    value <- .modifier_call_cdf(family, par, yy, ...)
    value[yy < lower] <- 0
    value[yy >= upper] <- 1
    if(!lower.tail) value <- 1 - value
    if(log.p) log(value) else value
  }
  fam$quantile <- function(par, p, lower.tail = TRUE, log.p = FALSE, ...) {
    probability <- if(log.p) exp(p) else p
    if(!lower.tail) probability <- 1 - probability
    pmin(pmax(family$quantile(par = par, p = probability, ...), lower), upper)
  }
  fam$random <- function(par, n, ...)
    pmin(pmax(family$random(par = par, n = n, ...), lower), upper)
  fam$mass <- function(par, y, log = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    value <- .modifier_mass(family, par, yy, log = TRUE, ...)
    if(is.finite(lower)) value[yy == lower] <- lower.logprob(par, n)[yy == lower]
    if(is.finite(upper)) value[yy == upper] <- upper.logprob(par, n)[yy == upper]
    value[yy < lower | yy > upper] <- -Inf
    if(log) value else exp(value)
  }
  fam$cdf_left <- function(par, y, log.p = FALSE, ...) {
    value <- fam$cdf(par, y, ...)
    value <- pmax(value - fam$mass(par, y, ...), 0)
    if(log.p) log(value) else value
  }

  links <- .modifier_links(family)
  fam$score <- list()
  for(parameter in family$names) {
    base.score <- family$score[[parameter]]
    fam$score[[parameter]] <- local({
      id <- parameter
      score <- base.score
      link <- links[[parameter]]
      function(par, y, ...) {
        n <- max(length(y), .modifier_rows(par))
        yy <- .modifier_recycle(y, n)
        value <- .modifier_recycle(score(par = par, y = yy, ...), n)
        if(is.finite(lower) && any(i <- yy == lower)) {
          pp <- .modifier_subset(par, i, n)
          fun <- function(p) lower.logprob(p, sum(i))
          value[i] <- .modifier_fd_score(fun, pp, id, link)
        }
        if(is.finite(upper) && any(i <- yy == upper)) {
          pp <- .modifier_subset(par, i, n)
          fun <- function(p) upper.logprob(p, sum(i))
          value[i] <- .modifier_fd_score(fun, pp, id, link)
        }
        value
      }
    })
  }
  fam$hessian <- .modifier_score_hessians(fam$score, family)
  fam$update <- NULL
  fam$log_likelihood <- function(par, y, ...)
    sum(fam$pdf(par = par, y = y, log = TRUE, ...), na.rm = TRUE)
  fam$valid.response <- function(y)
    all(is.na(y) | (y >= lower & y <= upper))
  fam$rqres <- function(par, y, ...) {
    lower.cdf <- fam$cdf_left(par, y)
    upper.cdf <- fam$cdf(par, y)
    stats::qnorm(stats::runif(length(y), lower.cdf, upper.cdf))
  }
  fam$support <- NULL
  fam$mean <- .modifier_no_moment("the censored mean")
  fam$variance <- .modifier_no_moment("the censored variance")
  class(fam) <- "gamlss2.family"
  fam
}

## Construct a family for observation-specific right censoring. The response
## is survival::Surv(time, status), with status equal to one for an event.
.right_censor_family <- function(family)
{
  response <- function(y) {
    if(!inherits(y, "Surv") || !identical(attr(y, "type"), "right"))
      stop("the response must be a right-censored 'Surv' object",
        call. = FALSE)
    list(time = as.numeric(y[, "time"]), status = as.numeric(y[, "status"]))
  }
  log_survival <- function(par, time)
    .modifier_call_cdf(
      family, par, time, lower.tail = FALSE, log.p = TRUE
    )

  fam <- .modifier_family(
    family, paste0("RightCensored(", family$family[1L], ")"), family$type
  )
  fam$pdf <- function(par, y, log = FALSE, ...) {
    yy <- response(y)
    value <- family$pdf(par = par, y = yy$time, log = TRUE, ...)
    censored <- yy$status == 0
    if(any(censored)) {
      n <- length(yy$time)
      value[censored] <- log_survival(
        .modifier_subset(par, censored, n), yy$time[censored]
      )
    }
    if(log) value else exp(value)
  }
  fam$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...)
    .modifier_call_cdf(
      family, par, y, lower.tail = lower.tail, log.p = log.p, ...
    )
  fam$quantile <- function(par, p, lower.tail = TRUE, log.p = FALSE, ...)
    .modifier_call_quantile(
      family, par, p, lower.tail = lower.tail, log.p = log.p, ...
    )
  fam$random <- function(par, n, ...)
    family$random(par = par, n = n, ...)

  links <- .modifier_links(family)
  fam$score <- list()
  for(parameter in family$names) {
    base.score <- family$score[[parameter]]
    fam$score[[parameter]] <- local({
      id <- parameter
      score <- base.score
      link <- links[[parameter]]
      function(par, y, ...) {
        yy <- response(y)
        n <- length(yy$time)
        value <- .modifier_recycle(
          score(par = par, y = yy$time, ...), n
        )
        censored <- yy$status == 0
        if(any(censored)) {
          pp <- .modifier_subset(par, censored, n)
          fun <- function(p)
            log_survival(p, yy$time[censored])
          value[censored] <- .modifier_fd_score(fun, pp, id, link)
        }
        value
      }
    })
  }
  fam$hessian <- .modifier_score_hessians(fam$score, family)
  fam$update <- NULL
  fam$log_likelihood <- function(par, y, ...)
    sum(fam$pdf(par = par, y = y, log = TRUE, ...), na.rm = TRUE)
  fam$valid.response <- function(y) {
    yy <- response(y)
    all(yy$status %in% 0:1) && family$valid.response(yy$time)
  }
  fam$rqres <- function(par, y, ...) {
    yy <- response(y)
    probability <- .modifier_call_cdf(family, par, yy$time, ...)
    censored <- yy$status == 0
    probability[censored] <- stats::runif(
      sum(censored), probability[censored], 1
    )
    stats::qnorm(probability)
  }
  if(!is.null(family$initialize)) {
    fam$initialize <- lapply(family$initialize, function(initialize) {
      if(!is.function(initialize)) return(initialize)
      function(y, ...) initialize(response(y)$time, ...)
    })
  }
  class(fam) <- "gamlss2.family"
  fam
}

## Add one or more explicit point masses to a base family.
inflate_family <- function(family = NO, at = 0)
{
  family <- complete_family(family)
  if(!is.numeric(at) || !length(at) || anyNA(at) ||
      any(!is.finite(at)) || anyDuplicated(at))
    stop("'at' must contain unique finite numeric values", call. = FALSE)
  at <- sort(as.numeric(at))
  k <- length(at)
  mixing.names <- paste0("pi", seq_len(k))
  if(any(mixing.names %in% family$names))
    stop("base and inflation parameter names overlap", call. = FALSE)

  parameter.rows <- function(par, n = NULL)
    .modifier_rows(as.list(par)[c(family$names, mixing.names)], n)
  probabilities <- function(par, n = NULL) {
    n <- parameter.rows(par, n)
    odds <- matrix(NA_real_, nrow = n, ncol = k)
    for(j in seq_len(k)) {
      value <- par[[mixing.names[j]]]
      if(is.null(value))
        stop("inflation parameter '", mixing.names[j], "' is missing",
          call. = FALSE)
      odds[, j] <- .modifier_recycle(value, n)
    }
    if(any(!is.finite(odds)) || any(odds < 0))
      stop("inflation odds must be finite and nonnegative", call. = FALSE)
    denominator <- 1 + rowSums(odds)
    cbind(base = 1 / denominator, odds / denominator)
  }
  base.parameters <- function(par)
    as.list(par)[family$names]
  atom.index <- function(y) match(y, at)

  fam <- .modifier_family(
    family,
    paste0("Inflated(", family$family[1L], ")"),
    if(identical(family$type, "continuous")) "mixed" else family$type
  )
  fam$names <- c(family$names, mixing.names)
  fam$links <- c(as.list(family$links),
    setNames(rep(list("log"), k), mixing.names))
  fam$map2par <- NULL

  fam$pdf <- function(par, y, log = FALSE, ...) {
    n <- parameter.rows(par, length(y))
    yy <- .modifier_recycle(y, n)
    pp <- base.parameters(par)
    probability <- probabilities(par, n)
    index <- atom.index(yy)
    value <- log(probability[, 1L]) +
      family$pdf(par = pp, y = yy, log = TRUE, ...)
    atoms <- !is.na(index)
    if(any(atoms)) {
      base.mass <- .modifier_mass(
        family, pp, yy, log = FALSE, ...
      )
      atom.probability <- probability[cbind(seq_len(n), index + 1L)]
      value[atoms] <- log(
        probability[atoms, 1L] * base.mass[atoms] +
          atom.probability[atoms]
      )
    }
    value <- as.numeric(value)
    value[is.na(yy)] <- NA_real_
    if(log) value else exp(value)
  }
  fam$mass <- function(par, y, log = FALSE, ...) {
    n <- parameter.rows(par, length(y))
    yy <- .modifier_recycle(y, n)
    pp <- base.parameters(par)
    probability <- probabilities(par, n)
    value <- probability[, 1L] *
      .modifier_mass(family, pp, yy, log = FALSE, ...)
    index <- atom.index(yy)
    atoms <- !is.na(index)
    if(any(atoms))
      value[atoms] <- value[atoms] +
        probability[cbind(which(atoms), index[atoms] + 1L)]
    value <- as.numeric(value)
    value[is.na(yy)] <- NA_real_
    if(log) log(value) else value
  }
  fam$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
    n <- parameter.rows(par, length(y))
    yy <- .modifier_recycle(y, n)
    probability <- probabilities(par, n)
    value <- probability[, 1L] * .modifier_call_cdf(
      family, base.parameters(par), yy, ...
    )
    for(j in seq_len(k))
      value <- value + probability[, j + 1L] * (yy >= at[j])
    value <- as.numeric(pmin(pmax(value, 0), 1))
    if(!lower.tail) value <- 1 - value
    if(log.p) log(value) else value
  }
  fam$cdf_left <- function(par, y, log.p = FALSE, ...) {
    value <- pmax(fam$cdf(par, y, ...) - fam$mass(par, y, ...), 0)
    if(log.p) log(value) else value
  }
  fam$quantile <- function(par, p, lower.tail = TRUE, log.p = FALSE, ...) {
    probability <- if(log.p) exp(p) else p
    if(!lower.tail) probability <- 1 - probability
    if(any(probability < 0 | probability > 1, na.rm = TRUE))
      stop("'p' must contain probabilities in [0, 1]", call. = FALSE)
    n <- parameter.rows(par, length(probability))
    probability <- .modifier_recycle(probability, n)
    component.probability <- probabilities(par, n)
    pp <- base.parameters(par)
    passed <- numeric(n)
    value <- rep.int(NA_real_, n)

    for(j in seq_len(k)) {
      left <- component.probability[, 1L] *
        .modifier_cdf_left(family, pp, rep.int(at[j], n)) + passed
      right <- component.probability[, 1L] *
        .modifier_call_cdf(family, pp, rep.int(at[j], n)) + passed +
        component.probability[, j + 1L]
      jump <- is.na(value) & probability > left & probability <= right
      value[jump] <- at[j]
      passed[probability > right] <- passed[probability > right] +
        component.probability[probability > right, j + 1L]
    }

    remaining <- is.na(value) & !is.na(probability)
    if(any(remaining)) {
      target <- (probability[remaining] - passed[remaining]) /
        component.probability[remaining, 1L]
      target <- pmin(pmax(target, 0), 1)
      value[remaining] <- family$quantile(
        par = .modifier_subset(pp, remaining, n), p = target, ...
      )
    }
    value
  }
  fam$random <- function(par, n, ...) {
    if(length(n) != 1L || !is.finite(n) || n < 0 || n != as.integer(n))
      stop("'n' must be a nonnegative integer", call. = FALSE)
    n <- as.integer(n)
    if(!n) return(numeric())
    probability <- probabilities(par, n)
    u <- stats::runif(n)
    component <- rowSums(u > t(apply(probability, 1L, cumsum)))
    value <- rep.int(NA_real_, n)
    base <- component == 0L
    if(any(base))
      value[base] <- family$random(
        par = .modifier_subset(base.parameters(par), base, n),
        n = sum(base), ...
      )
    for(j in seq_len(k))
      value[component == j] <- at[j]
    value
  }
  fam$probabilities <- function(par, n = NULL)
    as.data.frame(probabilities(par, n))
  if(is.function(family$mean)) {
    fam$mean <- function(par, ...) {
      n <- parameter.rows(par)
      probability <- probabilities(par, n)
      base.mean <- .modifier_recycle(family$mean(
        base.parameters(par), ...
      ), n)
      value <- probability[, 1L] * base.mean
      for(j in seq_len(k))
        value <- value + probability[, j + 1L] * at[j]
      as.numeric(value)
    }
  }
  if(is.function(family$mean) && is.function(family$variance)) {
    fam$variance <- function(par, ...) {
      n <- parameter.rows(par)
      probability <- probabilities(par, n)
      pp <- base.parameters(par)
      base.mean <- .modifier_recycle(family$mean(pp, ...), n)
      base.variance <- .modifier_recycle(family$variance(pp, ...), n)
      mixture.mean <- fam$mean(par, ...)
      second <- probability[, 1L] * (base.variance + base.mean^2)
      for(j in seq_len(k))
        second <- second + probability[, j + 1L] * at[j]^2
      as.numeric(pmax(second - mixture.mean^2, 0))
    }
  }

  base.initializers <- family$initialize
  fam$initialize <- list()
  for(parameter in family$names) {
    initialize <- base.initializers[[parameter]]
    if(is.function(initialize)) {
      fam$initialize[[parameter]] <- local({
        fun <- initialize
        function(y, ...) {
          keep <- !is.na(y) & !(y %in% at)
          if(any(keep)) fun(y[keep], ...) else fun(y, ...)
        }
      })
    }
  }
  for(j in seq_len(k)) {
    fam$initialize[[mixing.names[j]]] <- local({
      atom <- at[j]
      function(y, ...) {
        observed <- sum(y == atom, na.rm = TRUE)
        other <- sum(!is.na(y) & !(y %in% at))
        (observed + 0.5) / (other + 0.5)
      }
    })
  }

  responsibility <- function(par, y) {
    n <- parameter.rows(par, length(y))
    yy <- .modifier_recycle(y, n)
    probability <- probabilities(par, n)
    index <- atom.index(yy)
    base.responsibility <- rep.int(1, n)
    atom.responsibility <- matrix(0, nrow = n, ncol = k)
    atoms <- !is.na(index)
    if(any(atoms)) {
      base.mass <- .modifier_mass(
        family, base.parameters(par), yy, log = FALSE
      )
      atom.probability <- probability[cbind(seq_len(n), index + 1L)]
      total <- probability[, 1L] * base.mass + atom.probability
      base.responsibility[atoms] <-
        probability[atoms, 1L] * base.mass[atoms] / total[atoms]
      atom.responsibility[cbind(which(atoms), index[atoms])] <-
        atom.probability[atoms] / total[atoms]
    }
    list(base = base.responsibility, atom = atom.responsibility,
      probability = probability)
  }

  fam$score <- list()
  for(parameter in family$names) {
    base.score <- family$score[[parameter]]
    fam$score[[parameter]] <- local({
      score <- base.score
      function(par, y, ...) {
        n <- parameter.rows(par, length(y))
        yy <- .modifier_recycle(y, n)
        r <- responsibility(par, yy)$base
        value <- r * .modifier_recycle(score(
          par = base.parameters(par), y = yy, ...
        ), n)
        value[r == 0 & !is.finite(value)] <- 0
        value
      }
    })
  }
  for(j in seq_len(k)) {
    fam$score[[mixing.names[j]]] <- local({
      atom <- j
      function(par, y, ...) {
        r <- responsibility(par, y)
        r$atom[, atom] - r$probability[, atom + 1L]
      }
    })
  }
  fam$hessian <- .modifier_score_hessians(fam$score, fam)
  fam$update <- NULL
  fam$log_likelihood <- function(par, y, ...)
    sum(fam$pdf(par = par, y = y, log = TRUE, ...), na.rm = TRUE)
  fam$valid.response <- function(y) {
    keep <- is.na(y) | y %in% at
    if(all(keep)) return(TRUE)
    if(is.function(family$valid.response))
      family$valid.response(y[!keep])
    else TRUE
  }
  fam$rqres <- function(par, y, ...) {
    lower.cdf <- fam$cdf_left(par, y)
    upper.cdf <- fam$cdf(par, y)
    stats::qnorm(stats::runif(length(y), lower.cdf, upper.cdf))
  }
  fam$support <- NULL
  if(!is.function(fam$mean))
    fam$mean <- .modifier_no_moment("the inflated mean")
  if(!is.function(fam$variance))
    fam$variance <- .modifier_no_moment("the inflated variance")
  class(fam) <- "gamlss2.family"
  fam
}

## Add a hurdle and condition the base component away from that point.
hurdle_family <- function(family = PO, at = 0)
{
  family <- complete_family(family)
  if(!is.numeric(at) || length(at) != 1L || is.na(at) || !is.finite(at))
    stop("'at' must be one finite numeric value", call. = FALSE)

  log.normalizer <- function(par, n) {
    mass <- .modifier_mass(
      family, par, rep.int(at, n), log = FALSE
    )
    if(any(mass >= 1))
      stop("the base family is concentrated at the hurdle point",
        call. = FALSE)
    log1p(-mass)
  }
  adjustment <- function(par, y)
    -log.normalizer(par, max(length(y), .modifier_rows(par)))
  component <- .modifier_family(
    family, paste0("Excluded(", family$family[1L], ")")
  )
  component$pdf <- function(par, y, log = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    value <- family$pdf(par = par, y = yy, log = TRUE, ...) -
      log.normalizer(par, n)
    value[yy == at] <- -Inf
    if(log) value else exp(value)
  }
  component$mass <- function(par, y, log = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    value <- .modifier_mass(family, par, yy, log = TRUE, ...) -
      log.normalizer(par, n)
    value[yy == at] <- -Inf
    if(log) value else exp(value)
  }
  component$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    mass <- .modifier_mass(
      family, par, rep.int(at, n), log = FALSE
    )
    value <- (.modifier_call_cdf(family, par, yy, ...) -
      mass * (yy >= at)) / (1 - mass)
    value <- pmin(pmax(value, 0), 1)
    if(!lower.tail) value <- 1 - value
    if(log.p) log(value) else value
  }
  component$cdf_left <- function(par, y, log.p = FALSE, ...) {
    value <- pmax(component$cdf(par, y, ...) -
      component$mass(par, y, ...), 0)
    if(log.p) log(value) else value
  }
  component$quantile <- function(par, p, lower.tail = TRUE,
    log.p = FALSE, ...)
  {
    probability <- if(log.p) exp(p) else p
    if(!lower.tail) probability <- 1 - probability
    n <- max(length(probability), .modifier_rows(par))
    probability <- .modifier_recycle(probability, n)
    mass <- .modifier_mass(
      family, par, rep.int(at, n), log = FALSE
    )
    left <- .modifier_cdf_left(family, par, rep.int(at, n))
    split <- left / (1 - mass)
    target <- probability * (1 - mass)
    target[probability > split] <- target[probability > split] +
      mass[probability > split]
    family$quantile(par = par, p = pmin(pmax(target, 0), 1), ...)
  }
  component$random <- function(par, n, ...)
    component$quantile(par = par, p = stats::runif(n), ...)
  if(is.function(family$mean)) {
    component$mean <- function(par, ...) {
      n <- .modifier_rows(par)
      mass <- .modifier_mass(
        family, par, rep.int(at, n), log = FALSE
      )
      (family$mean(par, ...) - at * mass) / (1 - mass)
    }
  }
  if(is.function(family$mean) && is.function(family$variance)) {
    component$variance <- function(par, ...) {
      n <- .modifier_rows(par)
      mass <- .modifier_mass(
        family, par, rep.int(at, n), log = FALSE
      )
      base.mean <- .modifier_recycle(family$mean(par, ...), n)
      base.second <- .modifier_recycle(family$variance(par, ...), n) +
        base.mean^2
      conditional.mean <- component$mean(par, ...)
      pmax((base.second - at^2 * mass) / (1 - mass) -
        conditional.mean^2, 0)
    }
  }
  derivatives <- .modifier_add_derivatives(family, adjustment)
  component$score <- derivatives$score
  component$hessian <- derivatives$hessian
  component$update <- NULL
  component$initialize <- family$initialize
  component$log_likelihood <- function(par, y, ...)
    sum(component$pdf(par = par, y = y, log = TRUE, ...), na.rm = TRUE)
  component$rqres <- NULL
  component$support <- NULL
  if(!is.function(component$mean))
    component$mean <- .modifier_no_moment("the hurdle-component mean")
  if(!is.function(component$variance))
    component$variance <- .modifier_no_moment(
      "the hurdle-component variance"
    )
  class(component) <- "gamlss2.family"

  fam <- inflate_family(component, at = at)
  fam$family <- paste0("Hurdle(", family$family[1L], ")")
  fam$full.name <- "Hurdle distribution"
  fam
}

## Coarsen a base law into fixed adjacent intervals.
discretize_family <- function(family = NO, breaks,
  values = seq_len(length(breaks) - 1L), right = TRUE)
{
  family <- complete_family(family)
  if(!is.numeric(breaks) || length(breaks) < 2L || anyNA(breaks) ||
      is.unsorted(breaks, strictly = TRUE) ||
      breaks[1L] != -Inf || breaks[length(breaks)] != Inf) {
    stop("'breaks' must increase strictly from -Inf to Inf",
      call. = FALSE)
  }
  bins <- length(breaks) - 1L
  if(!is.numeric(values) || length(values) != bins || anyNA(values) ||
      any(!is.finite(values)) || is.unsorted(values, strictly = TRUE))
    stop("'values' must be finite, strictly increasing, and match the bins",
      call. = FALSE)
  if(!is.logical(right) || length(right) != 1L || is.na(right))
    stop("'right' must be TRUE or FALSE", call. = FALSE)

  bin.index <- function(y) match(y, values)
  interval.logprob <- function(par, index) {
    .modifier_log_interval(
      family, par, breaks[index], breaks[index + 1L],
      lower.left = !right, upper.left = !right
    )
  }
  fam <- .modifier_family(
    family, paste0("Discretized(", family$family[1L], ")"), "discrete"
  )
  fam$pdf <- function(par, y, log = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    index <- bin.index(yy)
    value <- rep.int(-Inf, n)
    ok <- !is.na(index)
    if(any(ok))
      value[ok] <- interval.logprob(
        .modifier_subset(par, ok, n), index[ok]
      )
    value[is.na(yy)] <- NA_real_
    if(log) value else exp(value)
  }
  fam$mass <- fam$pdf
  fam$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    index <- findInterval(yy, values)
    value <- numeric(n)
    inside <- index > 0L & index < bins
    if(any(inside)) {
      pp <- .modifier_subset(par, inside, n)
      boundary <- breaks[index[inside] + 1L]
      value[inside] <- if(right) {
        .modifier_call_cdf(family, pp, boundary)
      } else {
        .modifier_cdf_left(family, pp, boundary)
      }
    }
    value[index >= bins] <- 1
    value[is.na(yy)] <- NA_real_
    if(!lower.tail) value <- 1 - value
    if(log.p) log(value) else value
  }
  fam$cdf_left <- function(par, y, log.p = FALSE, ...) {
    value <- pmax(fam$cdf(par, y, ...) - fam$pdf(par, y, ...), 0)
    if(log.p) log(value) else value
  }
  fam$quantile <- function(par, p, lower.tail = TRUE, log.p = FALSE, ...) {
    probability <- if(log.p) exp(p) else p
    if(!lower.tail) probability <- 1 - probability
    if(any(probability < 0 | probability > 1, na.rm = TRUE))
      stop("'p' must contain probabilities in [0, 1]", call. = FALSE)
    n <- max(length(probability), .modifier_rows(par))
    probability <- .modifier_recycle(probability, n)
    value <- rep.int(values[bins], n)
    pending <- !is.na(probability)
    for(j in seq_len(bins)) {
      if(!any(pending)) break
      boundary <- rep.int(breaks[j + 1L], sum(pending))
      pp <- .modifier_subset(par, pending, n)
      distribution <- if(right) {
        .modifier_call_cdf(family, pp, boundary)
      } else {
        .modifier_cdf_left(family, pp, boundary)
      }
      take <- pending
      take[pending] <- probability[pending] <= distribution
      value[take] <- values[j]
      pending[take] <- FALSE
    }
    value[is.na(probability)] <- NA_real_
    value
  }
  fam$random <- function(par, n, ...)
    fam$quantile(par = par, p = stats::runif(n), ...)
  bin.probabilities <- function(par, n) {
    value <- matrix(NA_real_, nrow = n, ncol = bins)
    for(j in seq_len(bins))
      value[, j] <- exp(interval.logprob(par, rep.int(j, n)))
    value
  }
  fam$mean <- function(par, ...) {
    n <- .modifier_rows(par)
    as.numeric(bin.probabilities(par, n) %*% values)
  }
  fam$variance <- function(par, ...) {
    n <- .modifier_rows(par)
    probability <- bin.probabilities(par, n)
    center <- as.numeric(probability %*% values)
    as.numeric(probability %*% values^2 - center^2)
  }
  fam$score <- fam$hessian <- fam$update <- NULL
  fam$map2par <- family$map2par
  fam$log_likelihood <- function(par, y, ...)
    sum(fam$pdf(par = par, y = y, log = TRUE, ...), na.rm = TRUE)
  fam$valid.response <- function(y)
    all(is.na(y) | y %in% values)
  fam$rqres <- function(par, y, ...) {
    lower.cdf <- fam$cdf_left(par, y)
    upper.cdf <- fam$cdf(par, y)
    stats::qnorm(stats::runif(length(y), lower.cdf, upper.cdf))
  }
  fam$support <- NULL
  class(fam) <- "gamlss2.family"
  fam
}

## Round a continuous base law to a fixed number of decimal digits.
round_family <- function(family = NO, digits = 0L)
{
  family <- complete_family(family)
  if(!identical(tolower(family$type[1L]), "continuous"))
    stop("round_family() currently requires a continuous base family",
      call. = FALSE)
  if(!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits != as.integer(digits))
    stop("'digits' must be one integer", call. = FALSE)
  digits <- as.integer(digits)
  unit <- 10^(-digits)
  if(!is.finite(unit) || unit <= 0)
    stop("'digits' is outside the supported numeric range", call. = FALSE)
  on.grid <- function(y) {
    rounded <- round(y, digits)
    is.na(y) | abs(y - rounded) <=
      8 * .Machine$double.eps * pmax(1, abs(y))
  }

  fam <- .modifier_family(
    family, paste0("Rounded(", family$family[1L], ")"), "discrete"
  )
  fam$pdf <- function(par, y, log = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    value <- .modifier_log_interval(
      family, par, yy - unit / 2, yy + unit / 2
    )
    value[!on.grid(yy)] <- -Inf
    if(log) value else exp(value)
  }
  fam$mass <- fam$pdf
  fam$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
    n <- max(length(y), .modifier_rows(par))
    yy <- .modifier_recycle(y, n)
    grid <- floor(yy / unit + 8 * .Machine$double.eps) * unit
    value <- .modifier_call_cdf(family, par, grid + unit / 2)
    if(!lower.tail) value <- 1 - value
    if(log.p) log(value) else value
  }
  fam$cdf_left <- function(par, y, log.p = FALSE, ...) {
    value <- pmax(fam$cdf(par, y, ...) - fam$pdf(par, y, ...), 0)
    if(log.p) log(value) else value
  }
  fam$quantile <- function(par, p, lower.tail = TRUE, log.p = FALSE, ...) {
    probability <- if(log.p) exp(p) else p
    if(!lower.tail) probability <- 1 - probability
    round(family$quantile(par = par, p = probability, ...), digits)
  }
  fam$random <- function(par, n, ...)
    round(family$random(par = par, n = n, ...), digits)
  fam$score <- fam$hessian <- fam$update <- NULL
  fam$map2par <- family$map2par
  fam$log_likelihood <- function(par, y, ...)
    sum(fam$pdf(par = par, y = y, log = TRUE, ...), na.rm = TRUE)
  fam$valid.response <- function(y) all(on.grid(y))
  fam$rqres <- function(par, y, ...) {
    lower.cdf <- fam$cdf_left(par, y)
    upper.cdf <- fam$cdf(par, y)
    stats::qnorm(stats::runif(length(y), lower.cdf, upper.cdf))
  }
  fam$support <- NULL
  fam$mean <- .modifier_no_moment("the rounded mean")
  fam$variance <- .modifier_no_moment("the rounded variance")
  class(fam) <- "gamlss2.family"
  fam
}
