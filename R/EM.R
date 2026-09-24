## Expectation-maximization algorithm for mixture families.
EM <- function(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
{
  ## The optimizer is only available for families created by
  ## mixture_family().  Do not use the class alone here because all native
  ## gamlss2 families use the same class.
  is_mixture <- is.list(family$components) && length(family$components) > 1L &&
    is.list(family$component_map) &&
    is.function(family$probabilities) &&
    is.function(family$responsibilities)
  if(!is_mixture)
    stop("EM() can only be used with families returned by mixture_family()!")

  components <- family$components
  component_map <- family$component_map
  k <- length(components)
  n <- if(is.null(dim(y))) length(y) else nrow(y)

  ## Process offsets in the same form as RS().
  if(!is.null(offsets)) {
    if(is.null(nrow(offsets)) || nrow(offsets) < 1L)
      offsets <- NULL
    else
      offsets <- as.data.frame(offsets)
  }

  ## Parameter maps.
  component_parameters <- unlist(component_map, use.names = FALSE)
  if(!all(component_parameters %in% family$names))
    stop("invalid component parameter map in mixture family!")

  parameter_component <- rep(seq_len(k), lengths(component_map))
  parameter_inner <- unlist(lapply(component_map, names), use.names = FALSE)

  mixing_parameters <- setdiff(family$names, component_parameters)
  mixing_components <- suppressWarnings(as.integer(
    sub("^pi", "", mixing_parameters)
  ))
  if(length(mixing_parameters) != k - 1L ||
      anyNA(mixing_components) ||
      any(mixing_components < 1L | mixing_components > k) ||
      anyDuplicated(mixing_components)) {
    stop("invalid mixing parameter map in mixture family!")
  }
  names(mixing_components) <- mixing_parameters

  ## EM control parameters.
  eps <- control$em.eps
  if(is.null(eps))
    eps <- control$eps
  if(is.null(eps))
    eps <- 1e-05
  eps <- eps[1L]
  if(!is.numeric(eps) || length(eps) != 1L ||
      !is.finite(eps) || eps <= 0)
    stop("argument em.eps must be one positive finite value!")

  maxit <- control$em.maxit
  if(is.null(maxit))
    maxit <- control$maxit
  if(is.null(maxit))
    maxit <- 100L
  maxit <- maxit[1L]
  if(!is.numeric(maxit) || length(maxit) != 1L ||
      !is.finite(maxit) || maxit < 1L || maxit != floor(maxit))
    stop("argument em.maxit must be one positive integer!")
  maxit <- as.integer(maxit)

  mstep <- control$em.mstep
  if(is.null(mstep))
    mstep <- 1L
  if(!is.numeric(mstep) || length(mstep) != 1L ||
      !is.finite(mstep) || mstep < 1L || mstep != floor(mstep))
    stop("argument em.mstep must be one positive integer!")
  mstep <- as.integer(mstep)

  polish <- control$em.polish
  if(is.null(polish))
    polish <- TRUE
  polish <- isTRUE(polish)

  ## Links used for transforming component scores when a mixture parameter
  ## uses a link different from the corresponding component parameter.
  outer_links <- lapply(family$links, make.link2)
  component_links <- lapply(seq_len(k), function(component) {
    lapply(components[[component]]$links, make.link2)
  })

  normalize <- function(value, n) {
    value <- as.numeric(value)
    if(length(value) == n) value else rep(value, length.out = n)
  }

  component_par <- function(par, component, n) {
    map <- component_map[[component]]
    value <- lapply(map, function(parameter) {
      normalize(par[[parameter]], n)
    })
    names(value) <- names(map)
    value
  }

  ## Preserve coefficient starts and also provide parameter-level entries.
  ## The latter are needed by RS() when fixed parameters are present.
  warm_start <- function(fit) {
    value <- coef(fit, full = TRUE, lambdas = TRUE)
    missing <- setdiff(family$names, names(value))
    if(length(missing))
      value <- c(value, setNames(rep(NA_real_, length(missing)), missing))
    class(value) <- "coef.gamlss2"
    value
  }

  ## The observed mixture likelihood is used for checking convergence and for
  ## the final returned object.
  observed_loglik <- function(eta) {
    par <- family$map2par(eta)
    if(is.null(weights)) {
      family$log_likelihood(par = par, y = y)
    } else {
      sum(family$pdf(par = par, y = y, log = TRUE) * weights,
        na.rm = TRUE)
    }
  }

  ## Compute the observed-likelihood null deviance once.  The intermediate
  ## RS fits optimize the expected complete likelihood and therefore cannot
  ## supply this quantity.
  observed_null_deviance <- function()
  {
    eta <- initialize_eta(y, family, n, TRUE)
    beta <- vapply(family$names, function(parameter) {
      value <- eta[[parameter]]
      value <- value[is.finite(value)]
      if(length(value)) mean(value) else 0.0
    }, numeric(1L))

    fixed <- control$fixed
    if(is.null(fixed)) {
      fixed <- rep(FALSE, length(family$names))
      names(fixed) <- family$names
    } else {
      fixed <- as.list(fixed)
      if(is.null(names(fixed)))
        names(fixed) <- family$names[seq_along(fixed)]
      for(parameter in family$names) {
        if(is.null(fixed[[parameter]]))
          fixed[[parameter]] <- FALSE
      }
    }

    is_fixed <- function(parameter) {
      isTRUE(as.logical(fixed[[parameter]]))
    }

    for(parameter in family$names) {
      if(is_fixed(parameter)) {
        beta[parameter] <- make.link2(family$links[[parameter]])$linkfun(
          fixed[[parameter]]
        )
      }
    }

    null_eta <- function(value) {
      rval <- lapply(family$names, function(parameter) {
        rep(value[parameter], n)
      })
      names(rval) <- family$names
      if(!is.null(offsets)) {
        for(parameter in family$names) {
          if(!is.null(offsets[[parameter]]))
            rval[[parameter]] <- rval[[parameter]] + offsets[[parameter]]
        }
      }
      rval
    }

    ll <- observed_loglik(null_eta(beta))
    if(isTRUE(control$nullmodel) && length(unlist(xterms))) {
      objective <- function(value) {
        for(parameter in family$names) {
          if(is_fixed(parameter))
            value[parameter] <- beta[parameter]
        }
        value <- observed_loglik(null_eta(value))
        if(is.finite(value)) -value else .Machine$double.xmax
      }
      opt <- try(nlminb(beta, objective = objective), silent = TRUE)
      if(!inherits(opt, "try-error") && is.finite(opt$objective) &&
          -opt$objective > ll)
        ll <- -opt$objective
    }
    -2 * ll
  }

  ## Create the expected complete-data family for one EM iteration.  The
  ## responsibilities are initialized lazily from the predictors used by RS()
  ## and are then fixed for the complete M-step.
  em_family <- function()
  {
    state <- new.env(parent = emptyenv())
    state$responsibilities <- NULL
    state$logLik <- NULL

    responsibilities <- function(par, y) {
      if(is.null(state$responsibilities)) {
        state$logLik <- if(is.null(weights)) {
          family$log_likelihood(par = par, y = y)
        } else {
          sum(family$pdf(par = par, y = y, log = TRUE) * weights,
            na.rm = TRUE)
        }
        state$responsibilities <- as.matrix(
          family$responsibilities(par = par, y = y)
        )
        if(!identical(dim(state$responsibilities), c(n, k)))
          stop("invalid responsibility matrix returned by mixture family!")
      }
      state$responsibilities
    }

    ## Expected complete log-likelihood contribution for each observation.
    qpdf <- function(par, y, log = FALSE, ...)
    {
      posterior <- responsibilities(par, y)
      probabilities <- as.matrix(family$probabilities(par = par, n = n))
      if(!identical(dim(probabilities), c(n, k)))
        stop("invalid probability matrix returned by mixture family!")

      value <- matrix(0.0, nrow = n, ncol = k)
      for(component in seq_len(k)) {
        cp <- component_par(par, component, n)
        ld <- normalize(components[[component]]$pdf(
          par = cp, y = y, log = TRUE
        ), n)
        z <- posterior[, component] *
          (log(probabilities[, component]) + ld)
        z[posterior[, component] == 0 & !is.finite(z)] <- 0
        value[, component] <- z
      }
      value <- rowSums(value)
      if(log) value else exp(value)
    }

    qloglik <- function(par, y, ...) {
      sum(qpdf(par = par, y = y, log = TRUE, ...), na.rm = TRUE)
    }

    ## Working response and weights for the M-step. Component parameters use
    ## their ordinary GAMLSS working quantities multiplied by the fixed
    ## posterior probabilities. Mixing parameters use the multinomial-logit
    ## working quantities.
    qupdate <- function(par, y, eta, which)
    {
      posterior <- responsibilities(par, y)

      if(which %in% mixing_parameters) {
        component <- unname(mixing_components[which])
        probabilities <- as.matrix(family$probabilities(par = par, n = n))
        probability <- probabilities[, component]
        score <- posterior[, component] - probability
        hessian <- probability * (1 - probability)
      } else {
        occurrence <- which(component_parameters == which)
        if(!length(occurrence))
          stop("unknown mixture parameter in EM update: ", which, "!")

        outer_link <- outer_links[[which]]
        score <- hessian <- rep(0.0, n)

        ## A parameter can occur in more than one component. Summing the
        ## complete-data derivatives provides the shared-parameter M-step.
        for(index in occurrence) {
          component <- parameter_component[index]
          inner <- parameter_inner[index]
          cp <- component_par(par, component, n)
          score_fun <- components[[component]]$score[[inner]]
          hessian_fun <- components[[component]]$hessian[[inner]]
          if(!is.function(score_fun) || !is.function(hessian_fun))
            stop("component score or hessian is unavailable for parameter '",
              which, "'!")

          inner_link <- component_links[[component]][[inner]]
          same_link <- identical(outer_link$name, inner_link$name)

          if(same_link) {
            component_score <- normalize(score_fun(par = cp, y = y), n)
            component_hessian <- normalize(
              hessian_fun(par = cp, y = y), n
            )
            score <- score + posterior[, component] * component_score
            hessian <- hessian +
              posterior[, component] * component_hessian
          } else {
            ## This path is needed, for example, for the shifted links supplied
            ## by modal mixture initialization.
            score_at <- function(parameter_eta) {
              value <- outer_link$linkinv(parameter_eta)
              cpi <- cp
              cpi[[inner]] <- value
              component_eta <- inner_link$linkfun(value)
              multiplier <- outer_link$mu.eta(parameter_eta) /
                inner_link$mu.eta(component_eta)
              multiplier[!is.finite(multiplier)] <- 1
              posterior[, component] * normalize(
                score_fun(par = cpi, y = y), n
              ) * multiplier
            }

            parameter_eta <- outer_link$linkfun(par[[which]])
            step <- .Machine$double.eps^(1 / 3)
            component_score <- score_at(parameter_eta)
            score_plus <- score_at(parameter_eta + step)
            score_minus <- score_at(parameter_eta - step)
            score <- score + component_score
            hessian <- hessian -
              (score_plus - score_minus) / (2 * step)
          }
        }
      }

      score <- deriv_checks(score, is.weight = FALSE)
      hessian <- deriv_checks(hessian, is.weight = TRUE)
      list("eta" = eta + score / hessian, "weights" = hessian)
    }

    qfamily <- family
    qfamily$pdf <- qpdf
    qfamily$log_likelihood <- qloglik
    qfamily$update <- qupdate
    qfamily$.EM.state <- state
    qfamily
  }

  ## M-step controls. One complete RS sweep gives a generalized EM step.
  ## Additional sweeps can be requested with em.mstep.
  rs_control <- control
  rs_control$trace <- FALSE
  rs_control$flush <- FALSE
  rs_control$CG <- NULL
  rs_control$maxit <- c(mstep,
    if(is.null(control$maxit) || length(control$maxit) < 2L) {
      50L
    } else {
      control$maxit[2L]
    }
  )
  rs_control$nullmodel <- FALSE

  fit <- previous_fit <- NULL
  loglik <- numeric()
  converged <- FALSE
  stopped <- NULL
  em_start <- start

  for(iter in seq_len(maxit)) {
    qfamily <- em_family()
    candidate <- RS(x = x, y = y, specials = specials,
      family = qfamily, offsets = offsets, weights = weights,
      start = em_start, xterms = xterms, sterms = sterms,
      control = rs_control)

    candidate$family <- family
    ll <- observed_loglik(candidate$fitted.values)

    if(!is.finite(ll)) {
      if(is.null(previous_fit))
        stop("non-finite mixture log-likelihood in the first EM iteration!")
      stopped <- "non-finite mixture log-likelihood"
      fit <- previous_fit
      break
    }

    ## EM should not lower the observed likelihood. Retain the last accepted
    ## fit if numerical smoothing or an inexact M-step violates this property.
    if(length(loglik) && ll < tail(loglik, 1L)) {
      decrease <- tail(loglik, 1L) - ll
      tolerance <- eps * (abs(tail(loglik, 1L)) + 1e-08)
      if(decrease > tolerance) {
        stopped <- "mixture log-likelihood decreased"
        converged <- TRUE
        fit <- previous_fit
        break
      }
      stopped <- "numerical likelihood decrease"
      converged <- TRUE
      fit <- previous_fit
      break
    }

    previous_fit <- fit <- candidate
    loglik <- c(loglik, ll)

    initial_ll <- qfamily$.EM.state$logLik
    change <- abs(ll - initial_ll) / (abs(initial_ll) + 1e-08)
    if(!is.finite(change))
      change <- Inf

    if(control$trace) {
      cat(paste0("GAMLSS-EM iteration ", iter,
        ": Global Deviance = ", round(-2 * ll, digits = 4),
        " eps = ", fmt(change, width = 8, digits = 8), "    ",
        if(control$flush) "\r" else "\n"))
      if(control$flush)
        flush.console()
    }

    if(change <= eps) {
      converged <- TRUE
      break
    }

    ## Warm-start the next M-step, including smooth coefficients and smoothing
    ## parameters whenever they are available.
    em_start <- warm_start(fit)
  }

  if(control$trace && control$flush)
    cat("\n")

  if(is.null(fit))
    stop("EM algorithm did not produce a valid fit!")

  ## A final observed-likelihood RS fit is usually very short for parametric
  ## models.  It removes the linear-convergence tail of EM and, when accepted,
  ## supplies the ordinary fitted covariance information. Smooth models keep
  ## the EM fit because reselection of smoothing parameters can be expensive
  ## and need not increase the observed likelihood.
  polish_iterations <- 0L
  polish_loglik <- NA_real_
  polish_accepted <- FALSE
  polish_attempted <- polish && !length(unlist(sterms))
  if(polish_attempted) {
    polish_control <- control
    polish_control$trace <- FALSE
    polish_control$flush <- FALSE
    polish_control$nullmodel <- FALSE
    polish_control$eps <- eps
    polished <- RS(x = x, y = y, specials = specials,
      family = family, offsets = offsets, weights = weights,
      start = warm_start(fit),
      xterms = xterms, sterms = sterms, control = polish_control)
    polished$family <- family
    polished_loglik <- observed_loglik(polished$fitted.values)
    polish_iterations <- polished$iterations
    polish_loglik <- polished_loglik
    if(is.finite(polished_loglik) &&
        polished_loglik >= tail(loglik, 1L) - 1e-10 *
          (abs(tail(loglik, 1L)) + 1)) {
      fit <- polished
      polish_accepted <- TRUE
    }
  }

  if(identical(stopped, "non-finite mixture log-likelihood") &&
      !polish_accepted)
    warning(stopped, "; using the last accepted fit!")

  ## Restore observed-likelihood quantities.  The RS M-step reports the
  ## expected complete-data likelihood, which must not escape in the final
  ## gamlss2 object.
  final_loglik <- observed_loglik(fit$fitted.values)
  null_deviance <- observed_null_deviance()
  fit$logLik <- final_loglik
  fit$deviance <- -2 * final_loglik
  fit$null.deviance <- null_deviance
  fit$dev.reduction <- (null_deviance - fit$deviance) / null_deviance
  fit$nullmodel <- control$nullmodel
  fit$iterations <- length(loglik)
  fit$converged <- converged ||
    (polish_accepted && isTRUE(fit$converged))
  fit$control <- control
  fit$family <- family
  fit$posterior <- family$responsibilities(
    par = family$map2par(fit$fitted.values), y = y
  )
  fit$EM <- list(
    "logLik" = loglik,
    "iterations" = length(loglik),
    "converged" = converged,
    "mstep" = mstep,
    "polish" = polish,
    "polish.attempted" = polish_attempted,
    "polish.iterations" = polish_iterations,
    "polish.accepted" = polish_accepted,
    "polish.logLik" = polish_loglik,
    "stopped" = stopped
  )

  class(fit) <- "gamlss2"
  fit
}
