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

  if(!is.logical(initialize) || length(initialize) != 1L ||
      is.na(initialize)) {
    stop("'initialize' must be TRUE or FALSE", call. = FALSE)
  }

  components <- lapply(component_input, complete_family)
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

      initializer <- f$initialize[[inner]]
      if(initialize && inner %in% c("mu", "mean", "location")) {
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

  nonreference <- setdiff(seq_len(k), reference)
  mixing_names <- paste0("pi", nonreference)
  if(any(mixing_names %in% family_names))
    stop("component and mixing parameter names overlap", call. = FALSE)

  for(parameter in mixing_names) {
    links[[parameter]] <- "log"
    initializers[[parameter]] <- function(y, ...)
      rep(1, length(y))
  }
  family_names <- c(family_names, mixing_names)

  parameter_rows <- function(par, n = NULL) {
    lengths <- vapply(family_names, function(parameter) {
      value <- par[[parameter]]
      if(is.null(value)) 0L else length(value)
    }, integer(1L))
    max(c(1L, n, lengths))
  }

  component_parameters <- function(par, component, n) {
    map <- component_map[[component]]
    value <- lapply(unname(map), function(parameter) {
      z <- par[[parameter]]
      if(is.null(z))
        stop("parameter '", parameter, "' is missing", call. = FALSE)
      rep(z, length.out = n)
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
    yy <- rep(y, length.out = n)
    densities <- matrix(NA_real_, nrow = n, ncol = k)
    dots <- list(...)
    dots[c("par", "y", "log", "id")] <- NULL

    for(component in seq_len(k)) {
      component_par <- component_parameters(par, component, n)
      value <- do.call(
        components[[component]]$pdf,
        c(list(par = component_par, y = yy, log = TRUE), dots)
      )
      value <- rep(as.numeric(value), length.out = n)
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
    maxima <- apply(x, 1L, max)
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

  responsibility_matrix <- function(par, y, ...) {
    n <- parameter_rows(par, length(y))
    yy <- rep(y, length.out = n)
    joint <- log_joint_density(par, yy, ...)
    mixture <- row_log_sum_exp(joint, missing = is.na(yy))
    responsibilities <- exp(joint - mixture)

    invalid <- !is.finite(mixture) & !is.na(yy)
    if(any(invalid))
      responsibilities[invalid, ] <- probability_matrix(par, n)[invalid, ]

    colnames(responsibilities) <- component_labels
    responsibilities
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
          yy <- rep(y, length.out = n)
          responsibilities <- responsibility_matrix(par, yy, ...)
          component_par <- component_parameters(par, component_id, n)
          dots <- list(...)
          dots[c("par", "y", "id")] <- NULL
          score <- do.call(
            score_function,
            c(list(par = component_par, y = yy), dots)
          )
          score <- rep(as.numeric(score), length.out = n)
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
