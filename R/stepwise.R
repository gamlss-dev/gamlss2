## Helper functions.
get_df2 <- function(object) {
  df <- 0
  if(length(object$fitted.linear)) {
    if(is.null(attr(object$fitted.linear, "edf"))) {
      df <- df + sum(sapply(object$fitted.linear, function(x) length(x$coefficients)))
    } else {
      df <- df + sum(attr(object$fitted.linear, "edf"))
    }
  }
  if(length(object$fitted.specials)) {
    for(j in seq_along(object$fitted.specials)) {
      dfj <- sapply(object$fitted.specials[[j]], function(x) x$edf)
      df <- df + sum(unlist(dfj))
    }
  }
  return(df)
}

## Linear or nonlinear?
linear <- function(x) {
  grepl(".p.", x, fixed = TRUE)
}

## Helper function for splitting term names.
splitname <- function(x) {
  split <- if(linear(x)) ".p." else ".s."
  i <- regexpr(split, x, fixed = TRUE)[1L]
  return(c(substr(x, 1L, i - 1L),
    substr(x, i + nchar(split), nchar(x))))
}

## Helper function to add to xterms/sterms.
addterm <- function(x, y, term_map = NULL) {
  xi <- if(!is.null(term_map) && !is.null(term_map[[x]])) {
    term_map[[x]]
  } else {
    z <- splitname(x)
    list("parameter" = z[1L], "terms" = z[2L])
  }
  y[[xi$parameter]] <- unique(c(y[[xi$parameter]], xi$terms))
  return(y)
}

## Helper function to drop from xterms/sterms.
dropterm <- function(x, y, term_map = NULL) {
  xi <- if(!is.null(term_map) && !is.null(term_map[[x]])) {
    term_map[[x]]
  } else {
    z <- splitname(x)
    list("parameter" = z[1L], "terms" = z[2L])
  }
  y[[xi$parameter]] <- y[[xi$parameter]][!(y[[xi$parameter]] %in% xi$terms)]
  return(y)
}

## Extract formula terms separately for each distribution parameter.
formula_terms <- function(formula, parameters) {
  rval <- list(
    "linear" = setNames(vector("list", length(parameters)), parameters),
    "special" = setNames(vector("list", length(parameters)), parameters),
    "environment" = .GlobalEnv
  )
  for(i in parameters) {
    rval$linear[[i]] <- character(0)
    rval$special[[i]] <- character(0)
  }

  if(is.null(formula) || identical(formula, FALSE))
    return(rval)
  if(!inherits(formula, "formula") && !inherits(formula, "Formula") && !is.list(formula))
    stop("'fixed' must be NULL, a formula, a Formula object, or a list of formulas!")

  if(is.list(formula) && !inherits(formula, "Formula"))
    formula <- do.call("as.Formula", formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] > length(parameters))
    stop("'fixed' contains more predictors than the response distribution has parameters!")
  if(!is.null(environment(formula)))
    rval$environment <- environment(formula)

  for(j in seq_len(length(formula)[2L])) {
    fj <- formula(formula, lhs = 0, rhs = j)
    rhs <- paste(deparse(fj[[2L]]), collapse = "")
    if(identical(rhs, ".")) {
      rval$linear[[parameters[j]]] <- "."
      next
    }

    fl <- fake_formula(fj, nospecials = TRUE)
    rval$linear[[parameters[j]]] <- attr(terms(fl), "term.labels")
    rval$special[[parameters[j]]] <- fake_formula(fj, onlyspecials = TRUE)
  }
  return(rval)
}

## Split an interaction at top-level ':' operators.
interaction_parts <- function(x) {
  expr <- try(parse(text = x)[[1L]], silent = TRUE)
  if(inherits(expr, "try-error"))
    return(strsplit(x, ":", fixed = TRUE)[[1L]])

  split_expr <- function(z) {
    if(is.call(z) && identical(z[[1L]], as.name(":"))) {
      return(c(split_expr(z[[2L]]), split_expr(z[[3L]])))
    }
    paste(deparse(z), collapse = "")
  }
  return(split_expr(expr))
}

## Score the match between a formula term and a model-matrix column name.
term_match_score <- function(term, column) {
  term <- gsub("[ `]", "", term)
  column <- gsub("[ `]", "", column)
  tp <- interaction_parts(term)
  cp <- interaction_parts(column)
  if(length(tp) != length(cp))
    return(-Inf)

  score <- numeric(length(tp))
  for(j in seq_along(tp)) {
    if(identical(tp[j], cp[j])) {
      score[j] <- 100000L + nchar(tp[j])
    } else if(startsWith(cp[j], tp[j])) {
      score[j] <- nchar(tp[j])
    } else {
      return(-Inf)
    }
  }
  return(sum(score))
}

## Group model-matrix columns belonging to the same formula term. This is
## especially important for factors and interactions involving factors.
linear_term_groups <- function(labels, columns) {
  columns <- setdiff(columns, "(Intercept)")
  labels <- unique(setdiff(labels, c(".", "(Intercept)")))
  groups <- setNames(vector("list", length(labels)), labels)
  if(!length(columns))
    return(groups)

  if(length(labels)) {
    scores <- vapply(labels, function(label) {
      vapply(columns, function(column) term_match_score(label, column), numeric(1L))
    }, numeric(length(columns)))
    scores <- matrix(scores, nrow = length(columns), ncol = length(labels),
      dimnames = list(columns, labels))

    for(j in seq_along(columns)) {
      k <- which.max(scores[j, ])
      if(length(k) && is.finite(scores[j, k]))
        groups[[labels[k]]] <- c(groups[[labels[k]]], columns[j])
    }
  }

  used <- unique(unlist(groups, use.names = FALSE))
  for(column in setdiff(columns, used))
    groups[[column]] <- column
  groups[lengths(groups) > 0L]
}

## Build the selectable-term map and locate terms requested in fixed.
setup_stepwise_terms <- function(xterms, sterms, specials, parameters,
  scope = NULL, fixed = NULL)
{
  scope_terms <- formula_terms(scope, parameters)
  fixed_terms <- formula_terms(fixed, parameters)
  term_map <- list()

  for(parameter in parameters) {
    labels <- unique(c(scope_terms$linear[[parameter]], fixed_terms$linear[[parameter]]))
    groups <- linear_term_groups(labels, xterms[[parameter]])
    for(label in names(groups)) {
      id <- paste0(parameter, ".p.", label)
      term_map[[id]] <- list("parameter" = parameter, "type" = "p",
        "label" = label, "terms" = groups[[label]], "requires" = character(0))
    }

    if(length(sterms[[parameter]])) {
      slabels <- vapply(sterms[[parameter]], function(term) {
        z <- specials[[term]]
        if(is.list(z) && !is.null(z$orig.label)) z$orig.label else term
      }, character(1L))
      for(label in unique(slabels)) {
        id <- paste0(parameter, ".s.", label)
        term_map[[id]] <- list("parameter" = parameter, "type" = "s",
          "label" = label, "terms" = sterms[[parameter]][slabels == label],
          "requires" = character(0))
      }
    }
  }

  ## Linear interactions obey strong hierarchy whenever their marginal terms
  ## are present in the scope.
  normalize <- function(z) gsub("[ `]", "", z)
  for(id in names(term_map)) {
    z <- term_map[[id]]
    parts <- interaction_parts(z$label)
    if(z$type == "p" && length(parts) > 1L) {
      available <- names(term_map)[vapply(term_map, function(a) {
        identical(a$parameter, z$parameter) && identical(a$type, "p")
      }, logical(1L))]
      term_map[[id]]$requires <- available[vapply(term_map[available], function(a) {
        aparts <- interaction_parts(a$label)
        length(aparts) < length(parts) &&
          all(normalize(aparts) %in% normalize(parts))
      }, logical(1L))]
    }
  }

  fixed_ids <- character(0)
  for(parameter in parameters) {
    for(type in c("p", "s")) {
      requested <- fixed_terms[[if(type == "p") "linear" else "special"]][[parameter]]
      available <- names(term_map)[vapply(term_map, function(z) {
        identical(z$parameter, parameter) && identical(z$type, type)
      }, logical(1L))]
      if("." %in% requested) {
        fixed_ids <- c(fixed_ids, available)
        requested <- setdiff(requested, ".")
      }
      for(label in requested) {
        found <- available[vapply(term_map[available], function(z) {
          identical(gsub(" ", "", z$label), gsub(" ", "", label))
        }, logical(1L))]
        if(!length(found)) {
          stop(paste0("fixed term '", label,
            "' is not in the selection scope for parameter ", parameter, "!"))
        }
        fixed_ids <- c(fixed_ids, found)
      }
    }
  }

  ## Fixing an interaction also fixes its available marginal terms.
  repeat {
    required <- unique(unlist(lapply(term_map[fixed_ids], function(z) z$requires),
      use.names = FALSE))
    updated <- unique(c(fixed_ids, required))
    if(length(updated) == length(unique(fixed_ids)))
      break
    fixed_ids <- updated
  }

  return(list("term_map" = term_map, "fixed" = unique(fixed_ids),
    "environment" = if(is.null(fixed) || identical(fixed, FALSE)) {
      scope_terms$environment
    } else fixed_terms$environment))
}

term_parameter <- function(x, term_map) {
  if(!length(x)) return(character(0))
  vapply(term_map[x], function(z) z$parameter, character(1L))
}

term_label <- function(x, term_map = NULL) {
  if(!is.null(term_map) && !is.null(term_map[[x]])) {
    return(c(term_map[[x]]$parameter, term_map[[x]]$label))
  }
  splitname(x)
}

active_term_ids <- function(xterms, sterms, term_map) {
  if(is.null(term_map) || !length(term_map))
    return(character(0))
  ids <- names(term_map)
  ids[vapply(term_map, function(z) {
    current <- if(z$type == "p") xterms[[z$parameter]] else sterms[[z$parameter]]
    all(z$terms %in% current)
  }, logical(1L))]
}

valid_stepwise_terms <- function(xterms, sterms, term_map) {
  if(is.null(term_map))
    return(TRUE)
  active <- active_term_ids(xterms, sterms, term_map)
  all(vapply(term_map[active], function(z) all(z$requires %in% active), logical(1L)))
}

invalid_modelstats <- function() {
  list("df" = Inf, "GAIC" = Inf, "logLik" = -Inf)
}

## Helper function to extract GAIC etc.
modelstats <- function(model, K) {
  stats <- list("df" = get_df2(model))
  stats$GAIC <- model$deviance + K * stats$df
  stats$logLik <- model$logLik
  return(stats)
}

forward_backward_step <- function(selected, notselected, xterms, sterms, strategy, stats_save,
  x, y, specials, family, offsets, weights, start, control, trace, K,
  maxit = Inf, term_map = NULL, evaluate = NULL, remember = NULL)
{
  modelterms <- switch(strategy,
    "forward" = notselected,
    "backward" = selected
  )
  iter <- length(stats_save)
  cores <- control$cores

  if(is.null(evaluate)) {
    evaluate <- function(xtermsj, stermsj) {
      m <- RS(x, y, specials = specials, family, offsets,
        weights, start = start, xtermsj, sterms = stermsj, control)
      list("key" = NULL, "stats" = modelstats(m, K))
    }
  }
  if(is.null(remember))
    remember <- function(x) invisible(NULL)

  do <- TRUE
  k <- 0L
  while(do && (k < maxit)) {
    if(length(modelterms)) {
      fitfun <- function(j) {
        xtermsj <- xterms
        stermsj <- sterms

        if(strategy == "forward") {
          if(linear(j)) {
            xtermsj <- addterm(j, xtermsj, term_map)
          } else {
            stermsj <- addterm(j, stermsj, term_map)
          }
        } else {
          if(linear(j)) {
            xtermsj <- dropterm(j, xtermsj, term_map)
          } else {
            stermsj <- dropterm(j, stermsj, term_map)
          }
        }

        if(!valid_stepwise_terms(xtermsj, stermsj, term_map))
          return(list("key" = NULL, "stats" = invalid_modelstats()))
        evaluate(xtermsj, stermsj)
      }

      evaluated <- if(cores < 2L) {
        lapply(modelterms, fitfun)
      } else {
        parallel::mclapply(modelterms, fitfun, mc.cores = cores)
      }
      remember(evaluated)

      names(evaluated) <- modelterms
      stats <- lapply(evaluated, function(z) z$stats)
      gaic <- vapply(stats, function(z) z$GAIC, numeric(1L))
      j <- names(gaic)[which.min(gaic)]

      if(stats[[j]]$GAIC < stats_save[[iter]]$GAIC) {
        if(strategy == "forward") {
          if(linear(j)) {
            xterms <- addterm(j, xterms, term_map)
          } else {
            sterms <- addterm(j, sterms, term_map)
          }
          notselected <- notselected[notselected != j]
          selected <- c(selected, j)
        } else {
          if(linear(j)) {
            xterms <- dropterm(j, xterms, term_map)
          } else {
            sterms <- dropterm(j, sterms, term_map)
          }
          selected <- selected[selected != j]
          notselected <- c(notselected, j)
        }

        iter <- iter + 1L
        stats_save[[iter]] <- stats[[j]]
        stats_save[[iter]]$term <- j
        modelterms <- switch(strategy,
          "forward" = notselected,
          "backward" = selected
        )

        if(trace[2L]) {
          label <- term_label(j, term_map)
          cat("  GAIC = ", fmt(stats_save[[iter]]$GAIC, width = 10, digits = 4),
            if(strategy == "forward") " <+> parameter " else " <-> parameter ",
            label[1L], ", term ", label[2L],  "\n", sep = "")
        }
      } else {
        do <- FALSE
      }
    } else {
      do <- FALSE
    }
    k <- k + 1L
  }

  list("selected" = selected, "notselected" = notselected,
    "stats_save" = stats_save, "xterms" = xterms, "sterms" = sterms)
}

replace_step <- function(selected, notselected, xterms, sterms, stats_save,
  x, y, specials, family, offsets, weights, start, control, trace, K,
  maxit = Inf, term_map = NULL, evaluate = NULL, remember = NULL)
{
  iter <- length(stats_save)
  cores <- control$cores

  if(is.null(evaluate)) {
    evaluate <- function(xtermsj, stermsj) {
      m <- RS(x, y, specials = specials, family, offsets,
        weights, start = start, xtermsj, sterms = stermsj, control)
      list("key" = NULL, "stats" = modelstats(m, K))
    }
  }
  if(is.null(remember))
    remember <- function(x) invisible(NULL)

  do <- TRUE
  k <- 0L
  while(do && (k < maxit)) {
    if(length(selected) && length(notselected)) {
      pairs <- expand.grid("add" = notselected, "drop" = selected,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      fitfun <- function(k) {
        i <- pairs$add[k]
        j <- pairs$drop[k]
        xtermsj <- xterms
        stermsj <- sterms

        if(linear(j)) {
          xtermsj <- dropterm(j, xtermsj, term_map)
        } else {
          stermsj <- dropterm(j, stermsj, term_map)
        }
        if(linear(i)) {
          xtermsj <- addterm(i, xtermsj, term_map)
        } else {
          stermsj <- addterm(i, stermsj, term_map)
        }

        ev <- if(valid_stepwise_terms(xtermsj, stermsj, term_map)) {
          evaluate(xtermsj, stermsj)
        } else {
          list("key" = NULL, "stats" = invalid_modelstats())
        }
        ev$stats$drop <- j
        ev$stats$add <- i
        ev
      }

      evaluated <- if(cores < 2L) {
        lapply(seq_len(nrow(pairs)), fitfun)
      } else {
        parallel::mclapply(seq_len(nrow(pairs)), fitfun, mc.cores = cores)
      }
      remember(evaluated)

      stats <- lapply(evaluated, function(z) z$stats)
      best <- which.min(vapply(stats, function(z) z$GAIC, numeric(1L)))

      if(stats[[best]]$GAIC < stats_save[[iter]]$GAIC) {
        j <- stats[[best]]$drop
        i <- stats[[best]]$add

        if(linear(i)) {
          xterms <- addterm(i, xterms, term_map)
        } else {
          sterms <- addterm(i, sterms, term_map)
        }
        if(linear(j)) {
          xterms <- dropterm(j, xterms, term_map)
        } else {
          sterms <- dropterm(j, sterms, term_map)
        }

        selected <- selected[selected != j]
        notselected <- c(notselected, j)
        selected <- c(selected, i)
        notselected <- notselected[notselected != i]

        iter <- iter + 1L
        stats_save[[iter]] <- stats[[best]]
        stats_save[[iter]]$term <- i
        stats_save[[iter]]$dropped <- j

        if(trace[2L]) {
          dropped <- term_label(j, term_map)
          added <- term_label(i, term_map)
          cat("  GAIC = ", fmt(stats_save[[iter]]$GAIC, width = 10, digits = 4),
            " <-> parameter ", dropped[1L], ", term ", dropped[2L],  "\n", sep = "")
          cat("                    <+> parameter ", added[1L], ", term ",
            added[2L], "\n", sep = "")
        }
      } else {
        do <- FALSE
      }
    } else {
      do <- FALSE
    }
    k <- k + 1L
  }

  list("selected" = selected, "notselected" = notselected,
    "stats_save" = stats_save, "xterms" = xterms, "sterms" = sterms)
}

both_step <- function(selected, notselected, xterms, sterms, stats_save,
  x, y, specials, family, offsets, weights, start, control, trace, K,
  term_map = NULL, evaluate = NULL, remember = NULL)
{
  trace2 <- trace
  trace2[2L] <- FALSE
  gaic <- stats_save[[length(stats_save)]]$GAIC

  do <- TRUE
  while(do) {
    fs <- forward_backward_step(selected, notselected, xterms, sterms,
      strategy = "forward", stats_save, x, y, specials, family, offsets,
      weights, start, control, trace2, K, maxit = 1L, term_map = term_map,
      evaluate = evaluate, remember = remember)
    bs <- forward_backward_step(selected, notselected, xterms, sterms,
      strategy = "backward", stats_save, x, y, specials, family, offsets,
      weights, start, control, trace2, K, maxit = 1L, term_map = term_map,
      evaluate = evaluate, remember = remember)

    steps <- list("forward" = fs, "backward" = bs)
    step_gaic <- vapply(steps, function(z) {
      z$stats_save[[length(z$stats_save)]]$GAIC
    }, numeric(1L))
    best <- which.min(step_gaic)

    if(step_gaic[best] < gaic) {
      si <- steps[[best]]
      k <- length(si$stats_save)
      if(trace[2L]) {
        label <- term_label(si$stats_save[[k]]$term, term_map)
        cat("  GAIC = ", fmt(si$stats_save[[k]]$GAIC, width = 10, digits = 4),
          if(best == 1L) " <+> parameter " else " <-> parameter ",
          label[1L], ", term ", label[2L],  "\n", sep = "")
      }
      selected <- si$selected
      notselected <- si$notselected
      xterms <- si$xterms
      sterms <- si$sterms
      stats_save <- si$stats_save
      gaic <- step_gaic[best]
    } else {
      do <- FALSE
    }
  }

  list("selected" = selected, "notselected" = notselected,
    "stats_save" = stats_save, "xterms" = xterms, "sterms" = sterms)
}

all_step <- function(selected, notselected, xterms, sterms, stats_save,
  x, y, specials, family, offsets, weights, start, control, trace, K,
  term_map = NULL, evaluate = NULL, remember = NULL)
{
  trace2 <- trace
  trace2[2L] <- FALSE
  gaic <- stats_save[[length(stats_save)]]$GAIC

  do <- TRUE
  while(do) {
    fs <- forward_backward_step(selected, notselected, xterms, sterms,
      strategy = "forward", stats_save, x, y, specials, family, offsets,
      weights, start, control, trace2, K, maxit = 1L, term_map = term_map,
      evaluate = evaluate, remember = remember)
    bs <- forward_backward_step(selected, notselected, xterms, sterms,
      strategy = "backward", stats_save, x, y, specials, family, offsets,
      weights, start, control, trace2, K, maxit = 1L, term_map = term_map,
      evaluate = evaluate, remember = remember)
    rs <- replace_step(selected, notselected, xterms, sterms,
      stats_save, x, y, specials = specials, family, offsets, weights, start,
      control, trace2, K, maxit = 1L, term_map = term_map,
      evaluate = evaluate, remember = remember)

    steps <- list("forward" = fs, "backward" = bs, "replace" = rs)
    step_gaic <- vapply(steps, function(z) {
      z$stats_save[[length(z$stats_save)]]$GAIC
    }, numeric(1L))
    best <- which.min(step_gaic)

    if(step_gaic[best] < gaic) {
      si <- steps[[best]]
      k <- length(si$stats_save)
      if(trace[2L]) {
        if(best == 3L) {
          dropped <- term_label(si$stats_save[[k]]$dropped, term_map)
          added <- term_label(si$stats_save[[k]]$term, term_map)
          cat("  GAIC = ", fmt(si$stats_save[[k]]$GAIC, width = 10, digits = 4),
            " <-> parameter ", dropped[1L], ", term ", dropped[2L],  "\n", sep = "")
          cat("                    <+> parameter ", added[1L], ", term ",
            added[2L], "\n", sep = "")
        } else {
          label <- term_label(si$stats_save[[k]]$term, term_map)
          cat("  GAIC = ", fmt(si$stats_save[[k]]$GAIC, width = 10, digits = 4),
            if(best == 1L) " <+> parameter " else " <-> parameter ",
            label[1L], ", term ", label[2L],  "\n", sep = "")
        }
      }
      selected <- si$selected
      notselected <- si$notselected
      xterms <- si$xterms
      sterms <- si$sterms
      stats_save <- si$stats_save
      gaic <- step_gaic[best]
    } else {
      do <- FALSE
    }
  }

  list("selected" = selected, "notselected" = notselected,
    "stats_save" = stats_save, "xterms" = xterms, "sterms" = sterms)
}

stepwise <- function(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
{
  nx <- family$names

  ## Models must include intercepts.
  for(i in nx) {
    if(!("(Intercept)" %in% xterms[[i]]))
      stop(paste0("intercept is missing for parameter ", i, "!"))
  }

  ## Which strategies to use?
  choices <- c("forward.linear", "forward", "backward", "backward.linear",
    "replace", "replace.linear", "both", "both.linear", "all", "all.linear")
  strategy <- control$strategy
  if(is.null(strategy)) {
    strategy <- c("both.linear", "both")
  } else {
    strategy <- match.arg(strategy, choices, several.ok = TRUE)
  }

  keeporder <- isTRUE(control$keeporder)
  trace <- control$trace
  if(length(trace) < 2L)
    trace <- c(FALSE, trace)
  control$trace <- trace[1L]
  K <- if(is.null(control$K)) 2 else control$K
  if(is.null(control$cores))
    control$cores <- 1L

  ## A term can comprise multiple design-matrix columns, e.g., a factor or
  ## an interaction involving a factor.
  setup <- setup_stepwise_terms(xterms, sterms, specials, nx,
    scope = control$.stepwise_scope, fixed = control$.stepwise_fixed)
  term_map <- setup$term_map
  fixed <- setup$fixed
  candidates <- names(term_map)

  ## Begin with intercepts and the terms requested in fixed. Fixed terms are
  ## deliberately not put into selected, which is the set eligible for removal.
  xterms <- setNames(rep(list("(Intercept)"), length(nx)), nx)
  sterms <- setNames(rep(list(character(0)), length(nx)), nx)
  for(j in fixed) {
    if(linear(j)) {
      xterms <- addterm(j, xterms, term_map)
    } else {
      sterms <- addterm(j, sterms, term_map)
    }
  }
  selected <- character(0)
  notselected <- setdiff(candidates, fixed)

  ## Candidate models use a compact RS result and are cached by their exact
  ## ordered linear and special term sets. The ordered key preserves backfitting
  ## semantics while still recognizing revisited models.
  candidate_cache <- new.env(hash = TRUE, parent = emptyenv())
  candidate_control <- control
  candidate_control$.stepwise_candidate <- TRUE

  candidate_key <- function(xtermsj, stermsj) {
    parts <- unlist(lapply(nx, function(parameter) {
      c("parameter", parameter, "linear", xtermsj[[parameter]],
        "special", stermsj[[parameter]], "end")
    }), use.names = FALSE)
    paste(parts, collapse = "\034")
  }

  evaluate_candidate <- function(xtermsj, stermsj) {
    key <- candidate_key(xtermsj, stermsj)
    if(exists(key, envir = candidate_cache, inherits = FALSE))
      return(get(key, envir = candidate_cache, inherits = FALSE))

    m <- RS(x, y, specials = specials, family, offsets, weights,
      start = start, xtermsj, sterms = stermsj, candidate_control)
    value <- list("key" = key, "stats" = modelstats(m, K))
    assign(key, value, envir = candidate_cache)
    value
  }

  ## Forked workers update private copies of the cache, so merge their compact
  ## results into the parent cache after every candidate batch.
  remember_candidates <- function(values) {
    for(value in values) {
      if(is.null(value$key))
        next
      value$stats <- value$stats[c("df", "GAIC", "logLik")]
      assign(value$key, value, envir = candidate_cache)
    }
    invisible(NULL)
  }

  initial <- evaluate_candidate(xterms, sterms)
  remember_candidates(list(initial))
  stats_save <- list(initial$stats)

  ## Apply one search operation, optionally parameter by parameter.
  run_operation <- function(operation, candidate_ids, state) {
    run_subset <- function(ids, state) {
      isel <- intersect(state$selected, ids)
      inot <- intersect(state$notselected, ids)
      args <- list("evaluate" = evaluate_candidate, "remember" = remember_candidates)
      si <- switch(operation,
        "forward" = do.call(forward_backward_step, c(list(
          isel, inot, state$xterms, state$sterms, strategy = "forward",
          state$stats_save, x, y, specials, family, offsets, weights, start,
          control, trace, K, term_map = term_map), args)),
        "backward" = do.call(forward_backward_step, c(list(
          isel, inot, state$xterms, state$sterms, strategy = "backward",
          state$stats_save, x, y, specials, family, offsets, weights, start,
          control, trace, K, term_map = term_map), args)),
        "replace" = do.call(replace_step, c(list(
          isel, inot, state$xterms, state$sterms, state$stats_save,
          x, y, specials, family, offsets, weights, start, control, trace, K,
          term_map = term_map), args)),
        "both" = do.call(both_step, c(list(
          isel, inot, state$xterms, state$sterms, state$stats_save,
          x, y, specials, family, offsets, weights, start, control, trace, K,
          term_map = term_map), args)),
        "all" = do.call(all_step, c(list(
          isel, inot, state$xterms, state$sterms, state$stats_save,
          x, y, specials, family, offsets, weights, start, control, trace, K,
          term_map = term_map), args))
      )

      state$selected <- unique(c(setdiff(state$selected, ids), si$selected))
      state$notselected <- unique(c(setdiff(state$notselected, ids), si$notselected))
      state$xterms <- si$xterms
      state$sterms <- si$sterms
      state$stats_save <- si$stats_save
      state
    }

    if(keeporder) {
      for(parameter in if(operation == "backward") rev(nx) else nx) {
        ids <- candidate_ids[term_parameter(candidate_ids, term_map) == parameter]
        if(length(ids))
          state <- run_subset(ids, state)
      }
    } else if(length(candidate_ids)) {
      state <- run_subset(candidate_ids, state)
    }
    state
  }

  ## Individual operations already search to convergence. Multiple different
  ## operations still need complete sweeps because a later operation can open a
  ## move for an earlier one.
  run_phase <- function(operations, candidate_ids, state, suffix = "") {
    operation_order <- c("forward", "both", "all", "replace", "backward")
    operations <- operation_order[operation_order %in% operations]
    if(!length(operations))
      return(state)

    titles <- c("forward" = "Forward", "both" = "Bidirectional",
      "all" = "All Directions", "replace" = "Replace", "backward" = "Backward")
    sweep <- 1L
    repeat {
      gaic <- state$stats_save[[length(state$stats_save)]]$GAIC
      for(operation in operations) {
        if(trace[2L] && sweep == 1L)
          cat(titles[operation], suffix, " Selection\n", sep = "")
        state <- run_operation(operation, candidate_ids, state)
      }
      new_gaic <- state$stats_save[[length(state$stats_save)]]$GAIC
      if(length(operations) == 1L || !(new_gaic < gaic))
        break
      sweep <- sweep + 1L
      if(trace[2L])
        cat("Continue\n")
    }
    state
  }

  ## Move terms into the current model before a backward-only phase.
  initialize_upper <- function(ids, state, reset_stats = FALSE) {
    ids <- intersect(ids, state$notselected)
    for(j in ids) {
      if(linear(j)) {
        state$xterms <- addterm(j, state$xterms, term_map)
      } else {
        state$sterms <- addterm(j, state$sterms, term_map)
      }
    }
    state$selected <- unique(c(state$selected, ids))
    state$notselected <- setdiff(state$notselected, ids)

    upper <- evaluate_candidate(state$xterms, state$sterms)
    remember_candidates(list(upper))
    ms <- upper$stats
    if(reset_stats) {
      state$stats_save <- list(ms)
    } else {
      ms$term <- "(upper)"
      state$stats_save[[length(state$stats_save) + 1L]] <- ms
    }
    state
  }

  state <- list("selected" = selected, "notselected" = notselected,
    "xterms" = xterms, "sterms" = sterms, "stats_save" = stats_save)

  ## Preliminary search among linear terms only.
  linear_strategy <- sub(".linear", "", strategy[grepl(".linear", strategy, fixed = TRUE)],
    fixed = TRUE)
  linear_ids <- candidates[vapply(term_map, function(z) z$type == "p", logical(1L))]
  special_ids <- candidates[vapply(term_map, function(z) z$type == "s", logical(1L))]
  if("backward" %in% linear_strategy &&
      !any(c("forward", "both", "all") %in% linear_strategy)) {
    state <- initialize_upper(linear_ids, state, reset_stats = TRUE)
  }
  state <- run_phase(linear_strategy, linear_ids, state, suffix = " Linear")

  ## Search across linear and special terms. With a purely linear scope, do not
  ## repeat operations that the preliminary phase has already completed.
  full_strategy <- strategy[!grepl(".linear", strategy, fixed = TRUE)]
  if(!length(special_ids) && length(linear_strategy))
    full_strategy <- setdiff(full_strategy, linear_strategy)
  if("backward" %in% full_strategy &&
      !any(c("forward", "both", "all") %in% full_strategy)) {
    upper_ids <- if(length(linear_strategy)) special_ids else candidates
    state <- initialize_upper(upper_ids, state,
      reset_stats = length(state$stats_save) == 1L)
  }
  state <- run_phase(full_strategy, candidates, state)

  selected <- state$selected
  xterms <- state$xterms
  sterms <- state$sterms
  stats_save <- state$stats_save

  ## Always perform one complete final fit for inference and returned results.
  model <- RS(x, y, specials = specials, family, offsets,
    weights, start = start, xterms, sterms = sterms, control)

  ge <- function(j) {
    vapply(stats_save, function(z) as.numeric(z[[j]]), numeric(1L))
  }
  model$selection <- list("GAIC" = ge("GAIC"), "logLik" = ge("logLik"),
    "df" = ge("df"), "K" = K,
    "formula" = xs2formula(xterms, sterms, term_map, env = setup$environment))

  step_names <- vapply(seq_along(stats_save), function(j) {
    term <- stats_save[[j]]$term
    if(is.null(term))
      return(if(length(fixed)) "(fixed)" else "(Intercept)")
    if(!is.null(term_map[[term]])) {
      z <- term_map[[term]]
      return(paste(z$parameter, z$label, sep = "."))
    }
    term
  }, character(1L))
  names(model$selection$GAIC) <- step_names
  model$xterms <- xterms
  model$sterms <- sterms
  model$specials <- specials[unique(unlist(sterms))]

  if(trace[2L]) {
    cat("Final Model\n")
    for(j in names(model$selection$formula)) {
      cat(paste0("$", j, "\n.. "))
      print(model$selection$formula[[j]])
    }
  }

  return(model)
}

xs2formula <- function(x, s, term_map = NULL, env = .GlobalEnv)
{
  f <- list()
  parameters <- unique(c(names(x), names(s)))
  for(parameter in parameters) {
    if(is.null(term_map)) {
      sj <- s[[parameter]]
      if(length(sj) && any(grepl(").", sj, fixed = TRUE))) {
        sj <- vapply(strsplit(sj, ").", fixed = TRUE),
          function(z) paste0(z[1L], ")"), character(1L))
      }
      rhs <- c(unique(x[[parameter]]), unique(sj))
      rhs <- gsub("(Intercept)", "1", rhs, fixed = TRUE)
    } else {
      ids <- names(term_map)[vapply(term_map, function(z) {
        identical(z$parameter, parameter)
      }, logical(1L))]
      active <- ids[vapply(term_map[ids], function(z) {
        current <- if(z$type == "p") x[[parameter]] else s[[parameter]]
        all(z$terms %in% current)
      }, logical(1L))]
      labels <- unique(vapply(term_map[active], function(z) z$label, character(1L)))
      rhs <- c(if("(Intercept)" %in% x[[parameter]]) "1" else "0", labels)
    }
    if(!length(rhs))
      rhs <- "0"
    f[[parameter]] <- as.formula(paste("~", paste(rhs, collapse = "+")), env = env)
  }
  return(f)
}

step_gamlss2 <- function(formula, ..., fixed = NULL, K = 2,
  strategy = c("both.linear", "both"), keeporder = FALSE,
  cores = 1L)
{
  gamlss2(formula, ..., optimizer = stepwise, .stepwise_fixed = fixed,
    .stepwise_scope = formula, K = K, strategy = strategy,
    keeporder = keeporder, cores = cores)
}

new_formula <- function(object, ...)
{
  UseMethod("new_formula")
}

new_formula.default <- function(object, ...) {
  f <- object$selection$formula
  rn <- response_name(object)
  f[[1L]] <- eval(parse(text = paste("update(f[[1L]], ", rn, " ~ .)")))
  return(f)
}

