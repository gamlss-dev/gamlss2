## Asymmetric-Laplace quantile-regression family.
ALD <- function(tau = 0.5)
{
  if(length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1)
    stop("'tau' must be one finite number strictly between 0 and 1.")

  tau_info <- tau * (1 - tau)
  inv_tau_info <- 1 / tau_info

  ## Check / pinball loss.
  rho <- function(u) {
    u * (tau - (u < 0))
  }

  ## Subgradient wrt mu.
  psi_mu <- function(r) {
    tau - (r < 0) - 0.5 * (r == 0)
  }

  ## ALD quantile helper.
  qald <- function(p, mu, sigma) {
    n <- max(length(p), length(mu), length(sigma))

    p <- rep_len(p, n)
    mu <- rep_len(mu, n)
    sigma <- rep_len(sigma, n)

    if(any(!is.na(p) & (p < 0 | p > 1)))
      stop("probabilities must be in [0, 1].")

    out <- rep(NA_real_, n)

    lo <- !is.na(p) & p < tau
    hi <- !is.na(p) & !lo

    if(any(lo)) {
      out[lo] <- mu[lo] + sigma[lo] * (log(p[lo]) - log(tau)) / (1 - tau)
    }

    if(any(hi)) {
      out[hi] <- mu[hi] - sigma[hi] * (log1p(-p[hi]) - log1p(-tau)) / tau
    }

    out
  }

  fam <- list(
    family = "ALD",
    names = c("mu", "sigma"),
    links = c(
      mu = "identity",
      sigma = "log"
    ),
    type = "continuous",
    tau = tau,
    pdf = function(par, y, log = FALSE, ...) {
      u <- (y - par$mu) / par$sigma
      ans <-
        log(tau) +
        log1p(-tau) -
        log(par$sigma) -
        rho(u)

      if(log)
        ans
      else
        exp(ans)
    },
    log_likelihood = function(par, y, ...) {
      u <- (y - par$mu) / par$sigma
      ans <-
        log(tau) +
        log1p(-tau) -
        log(par$sigma) -
        rho(u)

      sum(ans, na.rm = TRUE)
    },
    cdf = function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
      u <- (y - par$mu) / par$sigma

      lp <- rep(NA_real_, length(u))

      lo <- !is.na(u) & u < 0
      hi <- !is.na(u) & !lo

      if(lower.tail) {
        if(any(lo)) {
          lp[lo] <-
            log(tau) +
            (1 - tau) * u[lo]
        }
        if(any(hi)) {
          lp[hi] <-
            log1p(
              -(1 - tau) *
              exp(-tau * u[hi])
            )
        }
      } else {
        if(any(lo)) {
          lp[lo] <-
            log1p(
              -tau *
              exp((1 - tau) * u[lo])
            )
        }
        if(any(hi)) {
          lp[hi] <-
            log1p(-tau) -
            tau * u[hi]
        }
      }
      if(log.p)
        lp
      else
        exp(lp)
    },
    quantile = function(par, p, lower.tail = TRUE, log.p = FALSE, ... ) {
      if(log.p) {
        p <- if(lower.tail)
          exp(p)
        else
          -expm1(p)
      } else if(!lower.tail) {
        p <- 1 - p
      }
      qald(
        p,
        par$mu,
        par$sigma
      )
    },
    random = function(n, par, ...) {
      n <- as.integer(n)
      qald(
        stats::runif(n),
        rep_len(par$mu, n),
        rep_len(par$sigma, n)
      )
    },
    mean = function(par, ...) {
      par$mu +
        par$sigma *
        (1 - 2 * tau) /
        (tau * (1 - tau))
    },
    variance = function(par, ...) {
      par$sigma^2 *
        (1 - 2 * tau + 2 * tau^2) /
        (tau^2 * (1 - tau)^2)
    },
    initialize = list(
      mu = function(y, ...) {
        q0 <- as.numeric(
          stats::quantile(
            y,
            probs = tau,
            type = 7,
            na.rm = TRUE,
            names = FALSE
          )
        )
        rep(q0, length(y))
      },
      sigma = function(y, ...) {
        q0 <- as.numeric(
          stats::quantile(
            y,
            probs = tau,
            type = 7,
            na.rm = TRUE,
            names = FALSE
          )
        )

        s0 <- mean(
          rho(y - q0),
          na.rm = TRUE
        )

        ys <- stats::IQR(
          y,
          na.rm = TRUE,
          type = 7
        )

        if(!is.finite(ys) || ys <= 0)
          ys <- stats::mad(
            y,
            constant = 1,
            na.rm = TRUE
          )

        if(!is.finite(ys) || ys <= 0)
          ys <- diff(
            range(y, na.rm = TRUE)
          )

        if(!is.finite(ys) || ys <= 0)
          ys <- max(1, abs(q0))

        tiny <- sqrt(.Machine$double.eps) * ys

        if(!is.finite(s0) || s0 <= tiny)
          s0 <- tiny

        rep(s0, length(y))
      }
    ),

    valid.response = function(y) {
      if(!is.numeric(y))
        stop("the response must be numeric.")
      if(all(is.na(y)))
        stop(
          "the response contains no non-missing values."
        )
      TRUE
    },

    score = list(
      mu = function(par, y, ...) {
        psi_mu(y - par$mu) / par$sigma
      },
      sigma = function(par, y, ...) {
        u <- (y - par$mu) / par$sigma
        rho(u) - 1
      }
    ),

    hessian = list(
      mu = function(par, y, ...) {
        tau_info / par$sigma^2
      },
      sigma = function(par, y, ...) {
        rep(1, length(y))
      }
    )
  )

  fam$update <- function(par, y, eta, which, ...) {
    if(identical(which, "mu")) {
      sigma <- par$sigma
      r <- y - par$mu
      psi <- tau - (r < 0) - 0.5 * (r == 0)
      return(list(
        eta = eta + sigma * psi * inv_tau_info,
        weights = tau_info / (sigma * sigma)
      ))
    }
    if(identical(which, "sigma")) {
      u <- (y - par$mu) / par$sigma
      return(list(
        eta = eta + u * (tau - (u < 0)) - 1,
        weights = rep.int(1, length(y))
      ))
    }
    stop("unknown ALD parameter: ", which)
  }

  class(fam) <- "gamlss2.family"

  fam
}

fit_ALD <- function(formula, data, tau,
  anchor = 0.5, max_step = 0.05,
  restarts = c("always", "on_failure", "never"),
  start = NULL, monotone = FALSE, verbose = interactive(), ...)
{
  cl <- match.call()
  restarts <- match.arg(restarts)

  if(!is.numeric(tau) || !length(tau) ||
      any(!is.finite(tau)) || any(tau <= 0 | tau >= 1))
    stop("'tau' must contain finite numbers strictly between 0 and 1.")

  if(length(anchor) != 1L || !is.finite(anchor) ||
      anchor <= 0 || anchor >= 1)
    stop("'anchor' must be one finite number strictly between 0 and 1.")

  if(length(max_step) != 1L || is.na(max_step) || max_step <= 0)
    stop("'max_step' must be one positive number.")

  if(length(monotone) != 1L || is.na(monotone))
    stop("'monotone' must be TRUE or FALSE.")
  monotone <- isTRUE(monotone)

  tau <- sort(unique(as.numeric(tau)))
  dots <- list(...)
  if("family" %in% names(dots))
    stop("'family' is supplied by fit_ALD().")

  ## Avoid printing the inner RS trace unless explicitly requested.
  if(is.null(dots$control) && is.null(dots$trace))
    dots$trace <- FALSE

  gamlss2_fun <- get("gamlss2", mode = "function", inherits = TRUE)
  tau_label <- function(x) {
    format(x, digits = 15L, trim = TRUE, scientific = FALSE)
  }

  fit_candidate <- function(q, candidate_start) {
    warns <- character()
    args <- c(
      list(
        formula = formula,
        data = data,
        family = ALD(tau = q),
        start = candidate_start
      ),
      dots
    )

    model <- withCallingHandlers(
      tryCatch(
        do.call(gamlss2_fun, args),
        error = function(e) e
      ),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    if(inherits(model, "error")) {
      return(list(
        valid = FALSE,
        model = NULL,
        logLik = -Inf,
        warnings = warns,
        error = conditionMessage(model)
      ))
    }

    ll <- tryCatch(as.numeric(stats::logLik(model)), error = function(e) NA_real_)
    pred <- tryCatch(
      stats::predict(model, model = "mu"),
      error = function(e) e
    )
    valid <- length(ll) == 1L && is.finite(ll) &&
      !inherits(pred, "error") && length(pred) && all(is.finite(pred))

    list(
      valid = valid,
      model = if(valid) model else NULL,
      logLik = if(valid) ll else -Inf,
      warnings = warns,
      error = if(valid) "" else "non-finite likelihood or predictions"
    )
  }

  fit_node <- function(q, previous = NULL, supplied_start = NULL,
    requested = FALSE)
  {
    candidate_start <- supplied_start
    start_name <- "start"

    if(!is.null(previous)) {
      candidate_start <- tryCatch(
        stats::coef(previous, full = TRUE, lambdas = TRUE),
        error = function(e) NULL
      )
      start_name <- "warm"
    }

    candidates <- list()
    if(!is.null(candidate_start))
      candidates[[start_name]] <- fit_candidate(q, candidate_start)

    need_cold <- is.null(candidate_start) ||
      identical(restarts, "always") ||
      (identical(restarts, "on_failure") &&
        !isTRUE(candidates[[start_name]]$valid))

    if(need_cold)
      candidates$cold <- fit_candidate(q, NULL)

    ok <- vapply(candidates, function(z) isTRUE(z$valid), logical(1L))
    if(!any(ok)) {
      detail <- vapply(
        candidates,
        function(z) z$error,
        character(1L)
      )
      stop(
        "all fits failed at tau = ", tau_label(q), ": ",
        paste(names(detail), detail, sep = ": ", collapse = "; ")
      )
    }

    valid_names <- names(candidates)[ok]
    chosen <- valid_names[
      which.max(vapply(candidates[valid_names], function(z) z$logLik, numeric(1L)))
    ]
    ans <- candidates[[chosen]]

    get_ll <- function(nm) {
      if(nm %in% names(candidates) && isTRUE(candidates[[nm]]$valid))
        candidates[[nm]]$logLik
      else
        NA_real_
    }

    diagnostic <- data.frame(
      tau = q,
      requested = requested,
      chosen = chosen,
      logLik = ans$logLik,
      warm.logLik = get_ll("warm"),
      start.logLik = get_ll("start"),
      cold.logLik = get_ll("cold"),
      iterations = as.integer(ans$model$iterations[1L]),
      warnings = paste(ans$warnings, collapse = " | "),
      stringsAsFactors = FALSE
    )

    if(isTRUE(verbose)) {
      message(
        "ALD tau = ", tau_label(q),
        " [", chosen, ", logLik = ",
        format(ans$logLik, digits = 8L), "]"
      )
    }

    list(model = ans$model, diagnostic = diagnostic)
  }

  is_requested <- function(q) {
    any(abs(tau - q) <= 1e-12 * pmax(1, abs(q)))
  }

  make_path <- function(targets, decreasing = FALSE) {
    if(!length(targets))
      return(numeric(0))

    targets <- sort(targets, decreasing = decreasing)
    current <- anchor
    path <- numeric(0)

    for(target in targets) {
      distance <- abs(target - current)
      nseg <- if(is.infinite(max_step)) {
        1L
      } else {
        max(1L, ceiling(distance / max_step - 1e-10))
      }
      nodes <- seq(current, target, length.out = nseg + 1L)[-1L]
      path <- c(path, nodes)
      current <- target
    }

    path
  }

  records <- list()
  anchor_fit <- fit_node(
    anchor,
    supplied_start = start,
    requested = is_requested(anchor)
  )
  records[[1L]] <- c(list(tau = anchor), anchor_fit)

  lower_path <- make_path(tau[tau < anchor], decreasing = TRUE)
  previous <- anchor_fit$model
  if(length(lower_path)) {
    for(q in lower_path) {
      node <- fit_node(q, previous = previous, requested = is_requested(q))
      records[[length(records) + 1L]] <- c(list(tau = q), node)
      previous <- node$model
    }
  }

  upper_path <- make_path(tau[tau > anchor], decreasing = FALSE)
  previous <- anchor_fit$model
  if(length(upper_path)) {
    for(q in upper_path) {
      node <- fit_node(q, previous = previous, requested = is_requested(q))
      records[[length(records) + 1L]] <- c(list(tau = q), node)
      previous <- node$model
    }
  }

  path_tau <- vapply(records, function(z) z$tau, numeric(1L))
  requested_index <- vapply(tau, function(q) {
    i <- which.min(abs(path_tau - q))
    if(abs(path_tau[i] - q) > 1e-10 * pmax(1, abs(q)))
      stop("internal error matching a requested quantile.")
    i
  }, integer(1L))

  models <- lapply(requested_index, function(i) records[[i]]$model)
  names(models) <- paste0("tau=", vapply(tau, tau_label, character(1L)))

  raw_fitted_values <- do.call(
    "cbind",
    lapply(models, stats::predict, model = "mu")
  )
  colnames(raw_fitted_values) <- names(models)
  fitted_values <- raw_fitted_values
  if(monotone && ncol(fitted_values) > 1L)
    fitted_values <- t(apply(fitted_values, 1L, sort))
  colnames(fitted_values) <- names(models)

  diagnostics <- do.call(
    "rbind",
    lapply(records, function(z) z$diagnostic)
  )
  rownames(diagnostics) <- NULL

  rval <- list(
    call = cl,
    tau = tau,
    models = models,
    monotone = monotone,
    raw.fitted.values = raw_fitted_values,
    fitted.values = fitted_values,
    diagnostics = diagnostics
  )
  class(rval) <- "ALDfit"
  rval
}


fitted.ALDfit <- function(object, ...)
{
  object$fitted.values
}


predict.ALDfit <- function(object, newdata = NULL, ...)
{
  dots <- list(...)
  dots[c("object", "model", "newdata")] <- NULL

  ans <- lapply(object$models, function(model) {
    args <- c(
      list(object = model, model = "mu"),
      if(is.null(newdata)) list() else list(newdata = newdata),
      dots
    )
    do.call(stats::predict, args)
  })

  ans <- do.call("cbind", ans)
  colnames(ans) <- names(object$models)
  if(isTRUE(object$monotone) && ncol(ans) > 1L)
    ans <- t(apply(ans, 1L, sort))
  colnames(ans) <- names(object$models)
  ans
}

print.ALDfit <- function(x, ...)
{
  requested <- x$diagnostics[x$diagnostics$requested, , drop = FALSE]
  requested <- requested[match(x$tau, requested$tau), , drop = FALSE]
  cat("ALD quantile fits:", length(x$tau), "\n")
  if(isTRUE(x$monotone))
    cat("Returned predictions use pointwise monotone rearrangement.\n")
  print(
    requested[, c("tau", "chosen", "logLik", "iterations", "warnings")],
    row.names = FALSE
  )
  invisible(x)
}

## if(FALSE) {
##   library("gamlss2")
##   data("mcycle", package = "MASS")

##   qu <- c(0.1, 0.5, 0.9)
##   m <- fit_ALD(
##     accel ~ s(times, k = 40, bs = "ad"), data = mcycle, tau = qu,
##     monotone = TRUE, verbose = FALSE
##   )
##   p <- fitted(m)

##   plot(accel ~ times, data = mcycle, ylim = range(mcycle$accel, p))
##   matlines(mcycle$times, p, col = 4, lty = 1, lwd = 2)
## }
