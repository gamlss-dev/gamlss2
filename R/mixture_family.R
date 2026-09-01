## Construct a finite mixture family.
mixture_family <- function(families = NO, k = NULL, reference = 1L,
  prefix = "c", initialize = TRUE)
{
  is_family_object <- function(x) {
    inherits(x, c("gamlss2.family", "gamlss.family", "distribution")) ||
      is.function(x) || is.character(x) ||
      (is.list(x) && (!is.null(x[["family"]]) || !is.null(x[["names"]]) ||
        !is.null(x[["pdf"]]) || !is.null(x[["d"]])))
  }

  is_component_list <- is.list(families) &&
    !is_family_object(families) &&
    length(families) > 0L

  if(is_component_list) {
    component_input <- families
    if(!is.null(k) && (!is.numeric(k) || length(k) != 1L ||
        !is.finite(k) || k != length(component_input))) {
      stop("'k' must equal the number of supplied component families",
        call. = FALSE)
    }
    k <- length(component_input)
  } else {
    if(is.null(k))
      k <- 2L
    if(!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
        k < 2 || k != as.integer(k)) {
      stop("'k' must be a single integer greater than or equal to 2",
        call. = FALSE)
    }
    k <- as.integer(k)
    component_input <- rep(list(families), k)
  }

  if(k < 2L)
    stop("at least two component families are required", call. = FALSE)

  if(!is.numeric(reference) || length(reference) != 1L ||
      !is.finite(reference) || reference != as.integer(reference) ||
      reference < 1L || reference > k) {
    stop("'reference' must identify one of the mixture components",
      call. = FALSE)
  }
  reference <- as.integer(reference)

  if(!is.character(prefix) || anyNA(prefix) ||
      !(length(prefix) %in% c(1L, k)) || any(!nzchar(prefix))) {
    stop("'prefix' must be one non-empty string or one per component",
      call. = FALSE)
  }
  component_labels <- if(length(prefix) == 1L) {
    paste0(prefix, seq_len(k))
  } else {
    prefix
  }
  if(anyDuplicated(component_labels))
    stop("component prefixes must be unique", call. = FALSE)

  if(is.logical(initialize) && length(initialize) == 1L &&
      !is.na(initialize)) {
    initialize_mode <- if(initialize) "separated" else "component"
  } else if(is.character(initialize) && length(initialize) == 1L &&
      !is.na(initialize)) {
    initialize_mode <- match.arg(
      initialize, c("separated", "cluster", "component")
    )
  } else {
    stop("'initialize' must be TRUE, FALSE, 'separated', 'cluster', or ",
      "'component'", call. = FALSE)
  }

  components <- if(is_component_list) {
    lapply(component_input, complete_family)
  } else {
    rep(list(complete_family(component_input[[1L]])), k)
  }
  if(any(!vapply(components, inherits, logical(1L),
      what = "gamlss2.family"))) {
    stop("all components must define gamlss2 families", call. = FALSE)
  }

  component_types <- vapply(components, function(x) {
    if(is.null(x$type)) "continuous" else tolower(x$type[1L])
  }, character(1L))
  if(length(unique(component_types)) != 1L) {
    stop("all mixture components must have the same distribution type",
      call. = FALSE)
  }

  component_map <- vector("list", k)
  family_names <- character()
  links <- list()
  initializers <- list()

  for(component in seq_len(k)) {
    f <- components[[component]]
    if(is.null(f$names) || !length(f$names))
      stop("each component family must have named parameters",
        call. = FALSE)

    inner_names <- as.character(f$names)
    outer_names <- paste0(component_labels[component], ".", inner_names)
    component_map[[component]] <- setNames(outer_names, inner_names)
    family_names <- c(family_names, outer_names)

    for(parameter in seq_along(inner_names)) {
      inner <- inner_names[parameter]
      outer <- outer_names[parameter]
      links[[outer]] <- f$links[[inner]]
    }
  }

  nonreference <- setdiff(seq_len(k), reference)
  mixing_names <- paste0("pi", nonreference)
  if(any(mixing_names %in% family_names))
    stop("component and mixing parameter names overlap", call. = FALSE)

  for(parameter in mixing_names)
    links[[parameter]] <- "log"
  family_names <- c(family_names, mixing_names)

  ## Deterministic response clustering for the optional "cluster" strategy.
  ## Midpoint quantiles define ordered centers and observations are assigned
  ## to the nearest center. A rank partition prevents empty clusters when
  ## quantiles coincide.
  response_partition <- function(y) {
    if(initialize_mode != "cluster" || !is.numeric(y) || !is.null(dim(y)))
      return(NULL)

    finite <- which(is.finite(y))
    if(!length(finite))
      return(NULL)
    yy <- as.numeric(y[finite])
    quantiles <- as.numeric(stats::quantile(
      yy, probs = (seq_len(k) - 0.5) / k,
      names = FALSE, type = 7
    ))

    centers <- quantiles
    cluster <- integer(length(yy))
    for(iteration in seq_len(25L)) {
      previous <- cluster
      distance <- abs(outer(yy, centers, "-"))
      cluster <- max.col(-distance, ties.method = "first")
      counts <- tabulate(cluster, nbins = k)

      if(length(yy) >= k && any(counts == 0L)) {
        ordering <- order(yy, seq_along(yy))
        ranked <- pmin(
          k,
          floor((seq_along(ordering) - 1L) * k / length(ordering)) + 1L
        )
        cluster[ordering] <- ranked
        counts <- tabulate(cluster, nbins = k)
      }
      for(component in seq_len(k)) {
        if(counts[component] > 0L)
          centers[component] <- mean(yy[cluster == component])
      }
      if(identical(cluster, previous) || any(counts == 0L))
        break
    }
    probabilities <- if(all(counts > 0L)) {
      counts / sum(counts)
    } else {
      rep(1 / k, k)
    }

    list(
      response = lapply(seq_len(k), function(component)
        y[finite[cluster == component]]),
      quantiles = quantiles,
      centers = centers,
      probabilities = probabilities
    )
  }

  initializer_result <- function(fun, y, ...) {
    if(!is.function(fun) || !length(y))
      return(NULL)
    value <- suppressWarnings(try(fun(y, ...), silent = TRUE))
    if(inherits(value, "try-error") && is.null(dim(y))) {
      value <- suppressWarnings(try(fun(matrix(y, ncol = 1L), ...),
        silent = TRUE))
    }
    if(inherits(value, "try-error"))
      return(NULL)
    value
  }

  ## Reduce a component-family initializer to one component-level value.
  initializer_value <- function(fun, y, ...) {
    value <- initializer_result(fun, y, ...)
    if(is.null(value))
      return(NULL)
    value <- suppressWarnings(as.numeric(value))
    value <- value[is.finite(value)]
    if(length(value)) mean(value) else NULL
  }

  ## Initial values are defined on the natural parameter scale. Only accept
  ## values that remain finite after applying the parameter link.
  valid_initial_values <- function(value, parameter) {
    value <- suppressWarnings(as.numeric(value))
    if(!length(value) || any(!is.finite(value)))
      return(FALSE)
    link <- make.link2(links[[parameter]])
    eta <- suppressWarnings(try(link$linkfun(value), silent = TRUE))
    !inherits(eta, "try-error") && length(eta) == length(value) &&
      all(is.finite(eta))
  }

  valid_initial_value <- function(value, parameter) {
    !is.null(value) && length(value) == 1L &&
      valid_initial_values(value, parameter)
  }

  ## If a cluster estimate lies on a link boundary, add one observation at
  ## the link-neutral value. This keeps log and logit starts in the interior.
  initial_candidates <- function(value, parameter, n = 0L) {
    if(is.null(value) || length(value) != 1L || !is.finite(value))
      return(list())
    candidates <- list(value)
    if(!valid_initial_value(value, parameter) && n > 0L) {
      link <- make.link2(links[[parameter]])
      neutral <- suppressWarnings(link$linkinv(0))
      regularized <- (n * value + neutral) / (n + 1)
      if(valid_initial_value(regularized, parameter))
        candidates[[2L]] <- regularized
    }
    candidates
  }

  component_initial_value <- function(y, component, inner, outer, ...) {
    partition <- response_partition(y)
    initializer <- components[[component]]$initialize[[inner]]
    candidates <- list()

    ## Preserve inherited observation-wise initialization for responses that
    ## cannot be partitioned, provided it is valid on the linked scale.
    if(is.null(partition)) {
      value <- initializer_result(initializer, y, ...)
      if(valid_initial_values(value, outer))
        return(rep(as.numeric(value), length.out = length(y)))
    }

    if(!is.null(partition)) {
      is_location <- inner %in% c("mu", "mean", "location")

      value <- initializer_value(
        initializer, partition$response[[component]], ...
      )
      candidates <- c(candidates, initial_candidates(
        value, outer, length(partition$response[[component]])
      ))
      if(is_location) {
        candidates <- c(candidates, initial_candidates(
          partition$centers[component], outer,
          length(partition$response[[component]])
        ))
        candidates <- c(candidates, initial_candidates(
          partition$quantiles[component], outer,
          length(partition$response[[component]])
        ))
      }
    }

    value <- initializer_value(initializer, y, ...)
    candidates <- c(candidates, initial_candidates(value, outer, length(y)))

    link <- make.link2(links[[outer]])
    candidates[[length(candidates) + 1L]] <-
      suppressWarnings(link$linkinv(0))
    for(value in candidates) {
      if(valid_initial_value(value, outer))
        return(rep(as.numeric(value), length(y)))
    }
    stop("could not find a finite initial value for parameter '", outer,
      "'", call. = FALSE)
  }

  for(component in seq_len(k)) {
    map <- component_map[[component]]
    for(inner in names(map)) {
      outer <- unname(map[[inner]])
      initializer <- components[[component]]$initialize[[inner]]
      if(initialize_mode == "cluster") {
        initializers[[outer]] <- local({
          component_id <- component
          inner_name <- inner
          outer_name <- outer
          function(y, ...) component_initial_value(
            y, component_id, inner_name, outer_name, ...
          )
        })
      } else if(initialize_mode == "separated" &&
          inner %in% c("mu", "mean", "location")) {
        probability <- (component - 0.5) / k
        initializers[[outer]] <- local({
          prob <- probability
          fallback <- initializer
          function(y, ...) {
            if(is.numeric(y)) {
              yy <- y[is.finite(y)]
              if(length(yy)) {
                value <- as.numeric(stats::quantile(
                  yy, probs = prob, names = FALSE, type = 7
                ))
                return(rep(value, length(y)))
              }
            }
            if(is.function(fallback))
              return(fallback(y, ...))
            rep(0, length(y))
          }
        })
      } else if(is.function(initializer)) {
        initializers[[outer]] <- initializer
      }
    }
  }

  for(j in seq_along(nonreference)) {
    component <- nonreference[j]
    parameter <- mixing_names[j]
    initializers[[parameter]] <- local({
      component_id <- component
      parameter_name <- parameter
      function(y, ...) {
        value <- 1
        if(initialize_mode == "cluster") {
          partition <- response_partition(y)
          if(!is.null(partition)) {
            value <- partition$probabilities[component_id] /
              partition$probabilities[reference]
          }
        }
        if(!valid_initial_value(value, parameter_name))
          value <- 1
        rep(value, length(y))
      }
    })
  }

  parameter_rows <- function(par, n = NULL) {
    max(c(1L, n, lengths(par[family_names])))
  }

  component_parameters <- function(par, component, n) {
    map <- component_map[[component]]
    value <- lapply(unname(map), function(parameter) {
      z <- par[[parameter]]
      if(is.null(z))
        stop("parameter '", parameter, "' is missing", call. = FALSE)
      if(length(z) == n) z else rep(z, length.out = n)
    })
    names(value) <- names(map)
    value
  }

  probability_matrix <- function(par, n = NULL) {
    n <- parameter_rows(par, n)
    odds <- matrix(1, nrow = n, ncol = k)

    for(j in seq_along(nonreference)) {
      component <- nonreference[j]
      value <- par[[mixing_names[j]]]
      if(is.null(value))
        stop("mixing parameter '", mixing_names[j], "' is missing",
          call. = FALSE)
      if(length(value) != n)
        value <- rep(value, length.out = n)
      if(any(!is.finite(value)) || any(value < 0))
        stop("mixing odds must be finite and nonnegative", call. = FALSE)
      odds[, component] <- value
    }

    probabilities <- odds / rowSums(odds)
    colnames(probabilities) <- paste0("pi", seq_len(k))
    probabilities
  }

  component_logdensities <- function(par, y, ...) {
    n <- parameter_rows(par, length(y))
    yy <- if(length(y) == n) y else rep(y, length.out = n)
    densities <- matrix(NA_real_, nrow = n, ncol = k)
    dots <- list(...)
    dots[c("par", "y", "log", "id")] <- NULL

    for(component in seq_len(k)) {
      component_par <- component_parameters(par, component, n)
      component_pdf <- components[[component]]$pdf
      value <- if(length(dots)) {
        do.call(
          component_pdf,
          c(list(par = component_par, y = yy, log = TRUE), dots)
        )
      } else {
        component_pdf(par = component_par, y = yy, log = TRUE)
      }
      value <- as.numeric(value)
      if(length(value) != n)
        value <- rep(value, length.out = n)
      invalid <- is.na(value) & !is.na(yy)
      value[invalid] <- -Inf
      densities[, component] <- value
    }

    colnames(densities) <- component_labels
    densities
  }

  log_joint_density <- function(par, y, ...) {
    n <- parameter_rows(par, length(y))
    probabilities <- probability_matrix(par, n)
    component_logdensities(par, y, ...) + log(probabilities)
  }

  row_log_sum_exp <- function(x, missing = NULL) {
    x[is.na(x)] <- -Inf
    maxima <- x[, 1L]
    if(ncol(x) > 1L) {
      for(column in 2:ncol(x))
        maxima <- pmax(maxima, x[, column])
    }
    value <- rep(-Inf, nrow(x))
    finite <- is.finite(maxima)
    if(any(finite)) {
      value[finite] <- maxima[finite] +
        log(rowSums(exp(x[finite, , drop = FALSE] - maxima[finite])))
    }
    if(!is.null(missing))
      value[missing] <- NA_real_
    value
  }

  responsibilities_from_parts <- function(logdensities, probabilities,
    missing = NULL)
  {
    joint <- logdensities + log(probabilities)
    mixture <- row_log_sum_exp(joint, missing = missing)
    responsibilities <- exp(joint - mixture)

    invalid <- !is.finite(mixture)
    if(!is.null(missing))
      invalid <- invalid & !missing
    if(any(invalid))
      responsibilities[invalid, ] <- probabilities[invalid, ]

    colnames(responsibilities) <- component_labels
    responsibilities
  }

  responsibility_matrix <- function(par, y, ...) {
    n <- parameter_rows(par, length(y))
    yy <- if(length(y) == n) y else rep(y, length.out = n)
    probabilities <- probability_matrix(par, n)
    logdensities <- component_logdensities(par, yy, ...)
    responsibilities_from_parts(
      logdensities, probabilities, missing = is.na(yy)
    )
  }

  component_names <- vapply(components, function(x) {
    as.character(x$family[1L])
  }, character(1L))

  fam <- list(
    family = paste0("Mixture(", paste(component_names, collapse = ", "), ")"),
    full.name = paste(k, "component finite mixture"),
    names = family_names,
    links = links,
    initialize = initializers,
    type = component_types[1L],
    components = components,
    component_map = component_map,
    reference = reference,
    mixing = setNames(nonreference, mixing_names)
  )

  fam$pdf <- function(par, y, log = FALSE, ...) {
    joint <- log_joint_density(par, y, ...)
    value <- row_log_sum_exp(
      joint,
      missing = is.na(rep(y, length.out = nrow(joint)))
    )
    if(log) value else exp(value)
  }

  fam$probabilities <- function(par, n = NULL, ...) {
    as.data.frame(probability_matrix(par, n))
  }

  fam$responsibilities <- function(par, y, ...) {
    as.data.frame(responsibility_matrix(par, y, ...))
  }

  fam$score <- list()
  for(component in seq_len(k)) {
    f <- components[[component]]
    map <- component_map[[component]]

    for(inner in names(map)) {
      outer <- unname(map[[inner]])
      component_score <- f$score[[inner]]
      fam$score[[outer]] <- local({
        component_id <- component
        score_function <- component_score
        function(par, y, ...) {
          n <- parameter_rows(par, length(y))
          yy <- if(length(y) == n) y else rep(y, length.out = n)
          responsibilities <- responsibility_matrix(par, yy, ...)
          component_par <- component_parameters(par, component_id, n)
          dots <- list(...)
          dots[c("par", "y", "id")] <- NULL
          score <- if(length(dots)) {
            do.call(
              score_function,
              c(list(par = component_par, y = yy), dots)
            )
          } else {
            score_function(par = component_par, y = yy)
          }
          score <- as.numeric(score)
          if(length(score) != n)
            score <- rep(score, length.out = n)
          value <- responsibilities[, component_id] * score
          value[responsibilities[, component_id] == 0 & !is.finite(score)] <- 0
          value
        }
      })
    }
  }

  for(j in seq_along(nonreference)) {
    component <- nonreference[j]
    parameter <- mixing_names[j]
    fam$score[[parameter]] <- local({
      component_id <- component
      function(par, y, ...) {
        n <- parameter_rows(par, length(y))
        responsibilities <- responsibility_matrix(par, y, ...)
        probabilities <- probability_matrix(par, n)
        responsibilities[, component_id] - probabilities[, component_id]
      }
    })
  }

  ## All parameter scores with a single responsibility calculation. This is
  ## especially useful for optimizing intercept-only mixture models, where a
  ## numerical gradient otherwise requires repeated full density evaluations.
  fam$score_matrix <- function(par, y, ...) {
    n <- parameter_rows(par, length(y))
    yy <- if(length(y) == n) y else rep(y, length.out = n)
    responsibilities <- responsibility_matrix(par, yy, ...)
    value <- matrix(NA_real_, nrow = n, ncol = length(family_names),
      dimnames = list(NULL, family_names))
    dots <- list(...)
    dots[c("par", "y", "id")] <- NULL

    for(component in seq_len(k)) {
      component_par <- component_parameters(par, component, n)
      map <- component_map[[component]]
      posterior <- responsibilities[, component]
      for(inner in names(map)) {
        outer <- unname(map[[inner]])
        score_function <- components[[component]]$score[[inner]]
        score <- if(length(dots)) {
          do.call(
            score_function,
            c(list(par = component_par, y = yy), dots)
          )
        } else {
          score_function(par = component_par, y = yy)
        }
        score <- as.numeric(score)
        if(length(score) != n)
          score <- rep(score, length.out = n)
        z <- posterior * score
        z[posterior == 0 & !is.finite(score)] <- 0
        value[, outer] <- z
      }
    }

    probabilities <- probability_matrix(par, n)
    for(j in seq_along(nonreference)) {
      component <- nonreference[j]
      value[, mixing_names[j]] <- responsibilities[, component] -
        probabilities[, component]
    }
    value
  }

  ## Compute mixture working responses and weights with a fused finite
  ## difference. The generic numerical Hessian evaluates all component
  ## densities for the center score and both perturbations. Here the center
  ## densities are shared, and only the changed component density is
  ## recomputed for a component-parameter perturbation. A mixing-parameter
  ## perturbation reuses all component densities. This is algebraically the
  ## same diagonal numerical Hessian used by complete_family().
  outer_parameters <- unlist(component_map, use.names = FALSE)
  parameter_component <- rep(seq_len(k), lengths(component_map))
  names(parameter_component) <- outer_parameters
  parameter_inner <- unlist(lapply(component_map, names), use.names = FALSE)
  names(parameter_inner) <- outer_parameters
  parameter_links <- lapply(links, make.link2)
  parameter_link_type <- vapply(links, function(link) {
    if(is.character(link) && length(link) == 1L &&
        link %in% c("identity", "log")) {
      link
    } else {
      "other"
    }
  }, character(1L))
  hessian_step <- .Machine$double.eps^(1 / 3)

  ## complete_family() otherwise builds a fully generic mapper that invokes a
  ## link closure for every parameter. Mixtures often contain many repeated
  ## identity and log links, so handle those common cases directly while
  ## retaining the exact generic sanitation for non-finite values.
  fam$map2par <- function(eta) {
    for(parameter in names(eta)) {
      value <- eta[[parameter]]
      link_type <- parameter_link_type[[parameter]]
      if(is.null(link_type))
        stop("unknown mixture parameter: ", parameter, call. = FALSE)
      if(link_type == "log") {
        value <- pmax(exp(value), .Machine$double.eps)
      } else if(link_type != "identity") {
        value <- parameter_links[[parameter]]$linkinv(value)
      }

      bad <- !is.finite(value)
      if(any(bad)) {
        index <- is.na(value)
        if(any(index))
          value[index] <- 0
        index <- value == Inf
        if(any(index))
          value[index] <- 10
        index <- value == -Inf
        if(any(index))
          value[index] <- -10
      }
      eta[[parameter]] <- value
    }
    eta
  }

  normalize_result <- function(value, n) {
    value <- as.numeric(value)
    if(length(value) == n) value else rep(value, length.out = n)
  }

  component_logdensity <- function(component, component_par, y) {
    value <- normalize_result(
      components[[component]]$pdf(
        par = component_par, y = y, log = TRUE
      ),
      length(y)
    )
    invalid <- is.na(value) & !is.na(y)
    value[invalid] <- -Inf
    value
  }

  weighted_component_score <- function(component, inner, component_par,
    y, posterior)
  {
    score <- normalize_result(
      components[[component]]$score[[inner]](
        par = component_par, y = y
      ),
      length(y)
    )
    value <- posterior * score
    value[posterior == 0 & !is.finite(score)] <- 0
    value
  }

  fam$update <- function(par, y, eta, which) {
    n <- parameter_rows(par, length(y))
    yy <- if(length(y) == n) y else rep(y, length.out = n)
    eta <- if(length(eta) == n) eta else rep(eta, length.out = n)
    missing <- is.na(yy)
    probabilities <- probability_matrix(par, n)
    logdensities <- component_logdensities(par, yy)
    responsibilities <- responsibilities_from_parts(
      logdensities, probabilities, missing
    )
    link <- parameter_links[[which]]
    if(is.null(link))
      stop("unknown parameter in update(): ", which, call. = FALSE)

    mixing_index <- match(which, mixing_names)
    if(!is.na(mixing_index)) {
      component <- nonreference[mixing_index]
      posterior <- responsibilities[, component]
      score <- posterior - probabilities[, component]

      parameter_eta <- link$linkfun(par[[which]])
      par_plus <- par_minus <- par
      par_plus[[which]] <- link$linkinv(parameter_eta + hessian_step)
      par_minus[[which]] <- link$linkinv(parameter_eta - hessian_step)

      probabilities_plus <- probability_matrix(par_plus, n)
      probabilities_minus <- probability_matrix(par_minus, n)
      responsibilities_plus <- responsibilities_from_parts(
        logdensities, probabilities_plus, missing
      )
      responsibilities_minus <- responsibilities_from_parts(
        logdensities, probabilities_minus, missing
      )
      score_plus <- responsibilities_plus[, component] -
        probabilities_plus[, component]
      score_minus <- responsibilities_minus[, component] -
        probabilities_minus[, component]
    } else {
      component <- unname(parameter_component[which])
      inner <- unname(parameter_inner[which])
      if(is.na(component) || is.na(inner))
        stop("unknown parameter in update(): ", which, call. = FALSE)

      component_par <- component_parameters(par, component, n)
      score <- weighted_component_score(
        component, inner, component_par, yy,
        responsibilities[, component]
      )

      parameter_eta <- link$linkfun(par[[which]])
      component_par_plus <- component_par_minus <- component_par
      component_par_plus[[inner]] <- rep(
        link$linkinv(parameter_eta + hessian_step), length.out = n
      )
      component_par_minus[[inner]] <- rep(
        link$linkinv(parameter_eta - hessian_step), length.out = n
      )

      logdensities_plus <- logdensities_minus <- logdensities
      logdensities_plus[, component] <- component_logdensity(
        component, component_par_plus, yy
      )
      logdensities_minus[, component] <- component_logdensity(
        component, component_par_minus, yy
      )
      responsibilities_plus <- responsibilities_from_parts(
        logdensities_plus, probabilities, missing
      )
      responsibilities_minus <- responsibilities_from_parts(
        logdensities_minus, probabilities, missing
      )
      score_plus <- weighted_component_score(
        component, inner, component_par_plus, yy,
        responsibilities_plus[, component]
      )
      score_minus <- weighted_component_score(
        component, inner, component_par_minus, yy,
        responsibilities_minus[, component]
      )
    }

    hessian <- -(score_plus - score_minus) / (2 * hessian_step)
    score <- deriv_checks(score, is.weight = FALSE)
    hessian <- deriv_checks(hessian, is.weight = TRUE)
    list(eta = eta + 1 / hessian * score, weights = hessian)
  }

  if(all(vapply(components, function(x) is.function(x$cdf),
      logical(1L)))) {
    fam$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
      n <- parameter_rows(par, length(y))
      yy <- rep(y, length.out = n)
      probabilities <- probability_matrix(par, n)
      values <- matrix(NA_real_, nrow = n, ncol = k)
      dots <- list(...)
      dots[c("par", "y", "lower.tail", "log.p", "log", "id")] <- NULL

      for(component in seq_len(k)) {
        component_par <- component_parameters(par, component, n)
        component_cdf <- components[[component]]$cdf
        cdf_args <- list(
          par = component_par,
          y = yy
        )
        cdf_formals <- names(formals(component_cdf))
        if("log.p" %in% cdf_formals)
          cdf_args$log.p <- FALSE
        else if("log" %in% cdf_formals)
          cdf_args$log <- FALSE
        values[, component] <- rep(as.numeric(do.call(
          component_cdf,
          c(cdf_args, dots)
        )), length.out = n)
      }

      value <- rowSums(probabilities * values)
      if(!lower.tail)
        value <- 1 - value
      if(log.p) log(value) else value
    }
  }

  if(is.function(fam$cdf) &&
      all(vapply(components, function(x) is.function(x$quantile),
        logical(1L)))) {
    fam$quantile <- function(par, p, lower.tail = TRUE,
      log.p = FALSE, ...)
    {
      n <- parameter_rows(par, length(p))
      probability <- rep(p, length.out = n)
      if(log.p)
        probability <- exp(probability)
      if(!lower.tail)
        probability <- 1 - probability
      if(any(probability < 0 | probability > 1, na.rm = TRUE))
        stop("'p' must contain probabilities in [0, 1]", call. = FALSE)

      missing <- is.na(probability)
      work_probability <- probability
      work_probability[missing] <- 0.5
      component_quantiles <- matrix(NA_real_, nrow = n, ncol = k)
      dots <- list(...)
      dots[c(
        "par", "p", "lower.tail", "log.p", "log", "id"
      )] <- NULL

      for(component in seq_len(k)) {
        component_par <- component_parameters(par, component, n)
        component_quantile <- components[[component]]$quantile
        quantile_args <- list(
          par = component_par,
          p = work_probability
        )
        quantile_formals <- names(formals(component_quantile))
        if("log.p" %in% quantile_formals)
          quantile_args$log.p <- FALSE
        else if("log" %in% quantile_formals)
          quantile_args$log <- FALSE
        component_quantiles[, component] <- rep(as.numeric(do.call(
          component_quantile,
          c(quantile_args, dots)
        )), length.out = n)
      }

      lower <- apply(component_quantiles, 1L, min)
      upper <- apply(component_quantiles, 1L, max)
      value <- rep(NA_real_, n)
      value[!missing & probability <= 0] <-
        lower[!missing & probability <= 0]
      value[!missing & probability >= 1] <-
        upper[!missing & probability >= 1]
      interior <- !missing & probability > 0 & probability < 1

      if(identical(component_types[1L], "discrete")) {
        lower <- floor(lower)
        upper <- ceiling(upper)
        for(iteration in seq_len(100L)) {
          active <- interior & lower < upper
          if(!any(active))
            break
          midpoint <- floor((lower + upper) / 2)
          distribution <- fam$cdf(par = par, y = midpoint)
          move_upper <- active & distribution >= probability
          move_lower <- active & !move_upper
          upper[move_upper] <- midpoint[move_upper]
          lower[move_lower] <- midpoint[move_lower] + 1
        }
        value[interior] <- lower[interior]
      } else {
        for(iteration in seq_len(80L)) {
          tolerance <- sqrt(.Machine$double.eps) *
            (1 + pmax(abs(lower), abs(upper)))
          active <- interior & is.finite(lower) & is.finite(upper) &
            (upper - lower) > tolerance
          if(!any(active))
            break
          midpoint <- (lower + upper) / 2
          distribution <- fam$cdf(par = par, y = midpoint)
          move_lower <- active & distribution < probability
          move_upper <- active & !move_lower
          lower[move_lower] <- midpoint[move_lower]
          upper[move_upper] <- midpoint[move_upper]
        }
        value[interior] <- (lower[interior] + upper[interior]) / 2
      }

      value
    }
  }

  if(all(vapply(components, function(x) is.function(x$random),
      logical(1L)))) {
    fam$random <- function(n, par, ...) {
      if(!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
          n < 0 || n != as.integer(n)) {
        stop("'n' must be a nonnegative integer", call. = FALSE)
      }
      n <- as.integer(n)
      if(n == 0L)
        return(numeric())

      probabilities <- probability_matrix(par, n)
      uniforms <- stats::runif(n)
      component_id <- 1L + rowSums(
        uniforms > t(apply(probabilities, 1L, cumsum))
      )
      value <- rep(NA_real_, n)
      dots <- list(...)
      dots[c("n", "par", "id")] <- NULL

      for(component in seq_len(k)) {
        index <- which(component_id == component)
        if(!length(index))
          next
        component_par <- component_parameters(par, component, n)
        component_par <- lapply(component_par, function(z) z[index])
        value[index] <- rep(as.numeric(do.call(
          components[[component]]$random,
          c(list(n = length(index), par = component_par), dots)
        )), length.out = length(index))
      }
      value
    }
  }

  if(all(vapply(components, function(x) is.function(x$mean),
      logical(1L)))) {
    fam$mean <- function(par, ...) {
      n <- parameter_rows(par)
      probabilities <- probability_matrix(par, n)
      values <- matrix(NA_real_, nrow = n, ncol = k)
      dots <- list(...)
      dots[c("par", "id")] <- NULL

      for(component in seq_len(k)) {
        component_par <- component_parameters(par, component, n)
        values[, component] <- rep(as.numeric(do.call(
          components[[component]]$mean,
          c(list(par = component_par), dots)
        )), length.out = n)
      }
      rowSums(probabilities * values)
    }
  }

  if(is.function(fam$mean) &&
      all(vapply(components, function(x) is.function(x$variance),
        logical(1L)))) {
    fam$variance <- function(par, ...) {
      n <- parameter_rows(par)
      probabilities <- probability_matrix(par, n)
      means <- variances <- matrix(NA_real_, nrow = n, ncol = k)
      dots <- list(...)
      dots[c("par", "id")] <- NULL

      for(component in seq_len(k)) {
        component_par <- component_parameters(par, component, n)
        means[, component] <- rep(as.numeric(do.call(
          components[[component]]$mean,
          c(list(par = component_par), dots)
        )), length.out = n)
        variances[, component] <- rep(as.numeric(do.call(
          components[[component]]$variance,
          c(list(par = component_par), dots)
        )), length.out = n)
      }

      mixture_mean <- rowSums(probabilities * means)
      pmax(
        rowSums(probabilities * (variances + means^2)) - mixture_mean^2,
        0
      )
    }
  }

  class(fam) <- "gamlss2.family"
  fam
}
