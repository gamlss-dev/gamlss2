cv_gamlss2 <- function(..., data, folds = 5,
  metric = log_pdf_metric, parallel = FALSE,
  simplify = TRUE)
{
  n <- nrow(data)

  if(!is.matrix(folds) & !is.list(folds) & !is.data.frame(folds)) {
    folds <- split(sample(seq_len(n)), rep_len(seq_len(folds), n))
  } else {
    if(!is.list(folds)) {
      folds <- as.list(as.data.frame(folds))
    }
  }

  runner <- function(i) {
    df <- data[-folds[[i]], , drop = FALSE]
    nd <- data[folds[[i]], , drop = FALSE]
    m <- gamlss2(..., data = df)
    list("score" = metric(model = m, data = nd), "id" = folds[[i]])
  }

  applier <- if(isTRUE(parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE))
      stop("parallel = TRUE requires the future.apply package.")
    future.apply::future_lapply
  } else {
    lapply
  }

  res <- applier(seq_along(folds), runner)

  if(simplify) {
    if(is.numeric(res[[1L]]$score) || is.vector(res[[1L]]$score)) {
      lens <- vapply(res, function(x) length(x$score), integer(1))
      if(all(lens < 2L)) {
        res <- vapply(res, function(x) as.numeric(x$score), numeric(1))
        names(res) <- paste0("fold_", seq_along(res))
      } else {
        res <- lapply(seq_along(res), function(i) {
          data.frame(id = res[[i]]$id, fold = i, score = res[[i]]$score)
        })
        res <- do.call("rbind", res)
        res <- res[order(res$id), c("fold", "score")]
        rownames(res) <- rownames(data)
      }
    }
  }

  return(res)
}

log_pdf_metric <- function(model, data) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::log_pdf(model, newdata = data)
}

rqres_metric <- function(model, data) {
  residuals(model, newdata = data)
}

mse_metric <- function(model, data) {
  ym <- mean(model, newdata = data)
  y <- data[[response_name(model)]]
  mean((y - ym)^2)
}

