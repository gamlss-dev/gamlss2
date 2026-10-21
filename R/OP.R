## Ordered probit family.
OP <- function(k) {
  if(length(k) != 1L || !is.numeric(k) || is.na(k) || !is.finite(k) ||
      k < 2 || k != floor(k) || k > .Machine$integer.max) {
    stop("argument k must be a single integer greater than or equal to 2.")
  }
  k <- as.integer(k)

  ## Parameter names: location and delta-encoded cutpoints.
  delta_indices <- if(k > 2L) seq.int(2L, k - 1L) else integer(0L)
  threshold_names <- c(
    "theta1", if(length(delta_indices)) paste0("delta", delta_indices)
  )
  par_names <- c("location", threshold_names)
  categories <- seq_len(k)
  category_levels <- as.character(categories)

  ## Identity links.
  links <- rep("identity", length(par_names))
  names(links) <- par_names

  as_categories <- function(y, argument = "y") {
    if(is.factor(y)) {
      if(!identical(levels(y), category_levels)) {
        stop(
          argument, " must have factor levels ",
          paste(category_levels, collapse = ", "), "."
        )
      }
      y <- as.integer(y)
    } else {
      if(!is.numeric(y))
        stop(argument, " must be numeric or a factor.")
      if(anyNA(y) || any(!is.finite(y)))
        stop(argument, " must not contain missing or non-finite values.")
      if(any(y != floor(y)))
        stop(argument, " must contain integer category values.")
      y <- as.integer(y)
    }

    if(any(y < 1L | y > k))
      stop(argument, " must contain category values in 1, ..., ", k, ".")

    y
  }

  compute_components <- function(par, n = NULL) {
    missing_parameters <- setdiff(par_names, names(par))
    if(length(missing_parameters)) {
      stop(
        "missing distribution parameter(s): ",
        paste(missing_parameters, collapse = ", "), "."
      )
    }

    parameter_lengths <- vapply(par[par_names], length, integer(1L))
    if(any(parameter_lengths < 1L))
      stop("distribution parameters must not be empty.")

    if(is.null(n))
      n <- max(parameter_lengths)
    else
      n <- max(c(n, parameter_lengths))

    if(any(!parameter_lengths %in% c(1L, n)))
      stop("distribution parameters must have length 1 or a common length.")

    ## Build increasing cutpoints.
    cuts <- matrix(NA_real_, nrow = n, ncol = k - 1L)
    cuts[, 1L] <- rep_len(par$theta1, n)
    if(k > 2L) {
      for(j in delta_indices) {
        cuts[, j] <- cuts[, j - 1L] +
          exp(rep_len(par[[paste0("delta", j)]], n))
      }
    }

    ## Cumulative probabilities: c_j = P(Y > j).
    cum_probs <- do.call(
      cbind,
      lapply(seq_len(k - 1L), function(j) {
        pnorm(rep_len(par$location, n) - cuts[, j])
      })
    )

    ## Category probabilities.
    probs <- matrix(NA_real_, nrow = n, ncol = k)
    probs[, 1L] <- 1 - cum_probs[, 1L]
    if(k > 2L) {
      for(j in delta_indices) {
        probs[, j] <- cum_probs[, j - 1L] - cum_probs[, j]
      }
    }
    probs[, k] <- cum_probs[, k - 1L]

    list(
      cuts = cuts,
      cum_probs = cum_probs,
      probs = probs
    )
  }

  initial_cutpoints <- function(y) {
    y <- as_categories(y, "response")
    if(!length(y))
      stop("response must not be empty.")

    counts <- tabulate(y, nbins = k)
    n <- sum(counts)
    cumulative_probs <- cumsum(counts)[seq_len(k - 1L)] / n

    ## Keep empirical cumulative probabilities away from 0 and 1 so
    ## initial cutpoints remain finite when boundary categories are empty.
    eps <- 0.5 / (n + 1)
    cumulative_probs <- pmin(pmax(cumulative_probs, eps), 1 - eps)
    qnorm(cumulative_probs)
  }

  ## Derivatives of the observed category probability on the predictor scale.
  compute_derivatives <- function(par, y, which, which2 = which) {
    comps <- compute_components(par, n = length(y))
    n <- nrow(comps$probs)
    if(length(y) == 1L && n > 1L)
      y <- rep(y, n)
    if(length(y) != n)
      stop("y must have length 1 or match the distribution parameters.")
    y_int <- as_categories(y)

    location <- rep_len(par$location, n)
    lower <- location - comps$cuts[cbind(seq_len(n), pmax(y_int - 1L, 1L))]
    upper <- location - comps$cuts[cbind(seq_len(n), pmin(y_int, k - 1L))]
    lower[y_int == 1L] <- Inf
    upper[y_int == k] <- -Inf
    f_lower <- dnorm(lower)
    f_upper <- dnorm(upper)
    g_lower <- -lower * f_lower
    g_upper <- -upper * f_upper
    g_lower[f_lower == 0] <- 0
    g_upper[f_upper == 0] <- 0

    ## Derivatives of location - gamma_j at the two category boundaries.
    slopes <- function(parameter) {
      if(parameter == "location")
        return(list(lower = 1, upper = 1))
      if(parameter == "theta1")
        return(list(lower = -1, upper = -1))
      j <- match(parameter, threshold_names)
      if(is.na(j))
        stop("unknown OP parameter: ", parameter)
      spacing <- exp(rep_len(par[[parameter]], n))
      dl <- du <- -spacing
      dl[y_int <= j | f_lower == 0] <- 0
      du[y_int < j | f_upper == 0] <- 0
      list(lower = dl, upper = du)
    }

    p <- comps$probs[cbind(seq_len(n), y_int)]
    active <- !is.na(p) & p >= 1e-8
    p <- pmax(p, 1e-8)
    d1 <- slopes(which)
    score <- (f_lower * d1$lower - f_upper * d1$upper) / p

    diagonal <- identical(which, which2)
    d2 <- if(diagonal) d1 else slopes(which2)
    score2 <- if(diagonal) score else
      (f_lower * d2$lower - f_upper * d2$upper) / p
    second <- (g_lower * d1$lower * d2$lower -
      g_upper * d1$upper * d2$upper) / p

    ## The second derivative of an exponential spacing equals its first.
    if(diagonal && which %in% threshold_names[-1L])
      second <- second + score
    hessian <- score * score2 - second

    ## The density is constant where its probability floor is active.
    score[!active] <- 0
    hessian[!active] <- 0
    list(score = score, hessian = hessian)
  }

  fam <- list(
    family = paste0("Ordered Probit (", k, " categories)"),
    names = par_names,
    links = links,

    pdf = function(par, y, log = FALSE, ...) {
      comps <- compute_components(par, n = length(y))
      n <- nrow(comps$probs)
      if(length(y) == 1L && n > 1L)
        y <- rep(y, n)
      if(length(y) != n)
        stop("y must have length 1 or match the distribution parameters.")
      y_int <- as_categories(y)

      p <- comps$probs[cbind(seq_len(n), y_int)]
      p[p < 1e-8 | is.na(p)] <- 1e-8

      if(log) {
        p <- log(p)
        p[is.na(p)] <- -1e10
      }

      p
    },

    initialize = {
      init_list <- list()

      ## Start with zero latent location so empirical cumulative probits
      ## directly initialize the cutpoints.
      init_list$location <- function(y, ...) {
        y <- as_categories(y, "response")
        rep(0, length(y))
      }

      ## Initialize theta1.
      init_list$theta1 <- function(y, ...) {
        cuts <- initial_cutpoints(y)
        rep(cuts[1L], length(y))
      }

      ## Initialize deltas as log-spacings between adjacent cutpoints.
      for(j in delta_indices) {
        init_list[[paste0("delta", j)]] <- local({
          jj <- j
          function(y, ...) {
            cuts <- initial_cutpoints(y)
            val <- log(max(cuts[jj] - cuts[jj - 1L], 1e-4))
            rep(val, length(y))
          }
        })
      }

      init_list
    }
  )

  ## Analytic scores and negative Hessians, including cross derivatives.
  fam$score <- fam$hessian <- list()
  for(i in seq_along(par_names)) {
    parameter <- par_names[i]
    fam$score[[parameter]] <- local({
      parameter <- par_names[i]
      function(par, y, ...) {
        compute_derivatives(par, y, parameter)$score
      }
    })
    fam$hessian[[parameter]] <- local({
      parameter <- par_names[i]
      function(par, y, ...) {
        compute_derivatives(par, y, parameter)$hessian
      }
    })
    if(i < length(par_names)) {
      for(j in seq.int(i + 1L, length(par_names))) {
        cross <- paste0(parameter, ":", par_names[j])
        fam$hessian[[cross]] <- local({
          parameter1 <- par_names[i]
          parameter2 <- par_names[j]
          function(par, y, ...) {
            compute_derivatives(par, y, parameter1, parameter2)$hessian
          }
        })
        fam$hessian[[paste0(par_names[j], ":", parameter)]] <- fam$hessian[[cross]]
      }
    }
  }

  ## Compute working response and weights with one probability calculation.
  fam$update <- function(par, y, eta, which, ...) {
    derivatives <- compute_derivatives(par, y, which)
    score <- deriv_checks(derivatives$score, is.weight = FALSE)
    hessian <- deriv_checks(derivatives$hessian, is.weight = TRUE)
    list(eta = eta + 1 / hessian * score, weights = hessian)
  }

  ## Probabilities on response scale.
  fam$probabilities <- function(par, ...) {
    probs <- compute_components(par)$probs
    colnames(probs) <- paste0("Pr(Y=", categories, ")")
    probs
  }

  ## Moments of the numeric category labels on the response scale.
  fam$mean <- function(par, ...) {
    probs <- fam$probabilities(par)
    drop(probs %*% categories)
  }

  fam$variance <- function(par, ...) {
    probs <- fam$probabilities(par)
    ey <- drop(probs %*% categories)
    drop(probs %*% categories^2) - ey^2
  }

  fam$transition <- function(par, ...) {
    tp <- compute_components(par)$cum_probs
    colnames(tp) <- paste0("Pr(Y>", seq_len(k - 1L), ")")
    tp
  }

  fam$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
    probs <- compute_components(par, n = length(y))$probs
    n <- nrow(probs)

    if(length(y) == 1L && n > 1L)
      y <- rep(y, n)
    if(length(y) != n)
      stop("y must have length 1 or match the distribution parameters.")
    y_int <- as_categories(y)

    cprobs <- t(apply(probs, 1L, cumsum))
    ans <- cprobs[cbind(seq_len(n), y_int)]

    if(!lower.tail)
      ans <- 1 - ans

    if(log.p)
      ans <- log(ans)

    ans
  }

  fam$quantile <- function(par, p, ...) {
    if(!is.numeric(p))
      stop("p must be numeric.")

    probs <- compute_components(par, n = length(p))$probs
    n <- nrow(probs)

    if(length(p) == 1L && n > 1L)
      p <- rep.int(p, n)

    if(length(p) != n)
      stop("p must have length 1 or match the distribution parameters.")
    if(anyNA(p) || any(!is.finite(p)))
      stop("p must not contain missing or non-finite values.")
    if(any(p < 0 | p > 1))
      stop("p must be in [0, 1].")

    cprobs <- t(apply(probs, 1L, cumsum))

    q <- integer(n)
    for(i in seq_len(n)) {
      idx <- which(cprobs[i, ] >= p[i])[1L]
      if(is.na(idx))
        idx <- k
      q[i] <- idx
    }

    q
  }

  fam$log_likelihood <- function(par, y, ...) {
    sum(fam$pdf(par, y, log = TRUE))
  }

  fam$valid.response <- function(x) {
    as_categories(x, "response")
    TRUE
  }

  fam$residuals <- function(object, ...) {
    rqres_oprobit(object, ...)
  }

  fam$type <- "discrete"

  class(fam) <- c("gamlss2.family", "family.bamlss")
  fam
}

rqres_oprobit <- function(object, ...) {
  fam <- family(object)
  if(!grepl("^Ordered Probit", fam$family))
    stop("OP family required.")

  mf <- model.frame(object)
  y <- stats::model.response(mf)

  par <- predict(object, type = "parameter")
  probs <- fam$probabilities(par)

  K <- ncol(probs)
  y_int <- if(is.factor(y)) {
    match(as.character(y), as.character(seq_len(K)))
  } else {
    as.integer(y)
  }
  n <- length(y_int)

  if(any(y_int < 1L | y_int > K | is.na(y_int))) {
    bad_vals <- sort(unique(y[is.na(y_int) | y_int < 1L | y_int > K]))
    stop(
      "response has values outside 1..", K, " implied by OP(k).\n",
      "Offending values: ", paste(bad_vals, collapse = ", ")
    )
  }

  cprobs <- t(apply(probs, 1L, cumsum))
  F_upper <- cprobs[cbind(seq_len(n), y_int)]

  F_lower <- numeric(n)
  idx_gt1 <- which(y_int > 1L)
  if(length(idx_gt1) > 0L) {
    F_lower[idx_gt1] <- cprobs[cbind(idx_gt1, y_int[idx_gt1] - 1L)]
  }

  u <- stats::runif(n, min = F_lower, max = F_upper)
  u[u <= 0] <- 1e-12
  u[u >= 1] <- 1 - 1e-12
  stats::qnorm(u)
}

oprobit <- function(k) {
  OP(k = k)
}
