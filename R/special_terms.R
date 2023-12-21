## Special conditional inference forest constructor.
cf <- function(formula, ...)
{
  stopifnot(requireNamespace("partykit"))
  st <- list()
  ctr <- list(...)
  ntree <- ctr$ntree
  if(is.null(ntree))
    ntree <- 100L
  ctr$ntree <- NULL
  st$control <- do.call(partykit::ctree_control, ctr)
  st$formula <- formula
  st$term <- all.vars(formula)
  st$label <- paste0("cf(", paste0(gsub(" ", "", as.character(formula)), collapse = ""), ")")
  st$data <- model.frame(formula)
  st$ntree <- ntree
  class(st) <- c("special", "cf")
  return(st)
}

## The fitting function for the backfitting algorithm.
special_fit.cf <- function(x, z, w, y, eta, j, family, control, ...)
{
  f <- update(x$formula, response_z ~ .)
  x$data$response_z <- z
  rval <- list(
    "model" = partykit::cforest(formula = f, data = x$data, weights = w,
      control = x$control, ntree = x$ntree)
  )
  rval$fitted.values <- predict(rval$model)
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift
  rval$edf <- x$ntree
  class(rval) <- "cf.fitted"
  return(rval)
}

## A predict method.
special_predict.cf.fitted <- function(x, data, se.fit = FALSE, ...)
{
  if(se.fit) {
    ct <- lapply(seq_along(x$model$nodes), function(i) partykit::as.constparty(
      partykit::party(x$model$nodes[[i]], data = x$model$data, terms = x$model$terms,
      fitted = data.frame(
        `(response)` = x$model$fitted[["(response)"]],
        `(weights)` = x$model$weights[[i]],
        check.names = FALSE))
    ))
    p <- sapply(ct, predict, newdata = data)
    p <- apply(p, 1, quantile, prob = c(0.05, 0.5, 0.95))
    p <- t(p)
    colnames(p) <- c("lower", "fit", "upper")
  } else {
    p <- predict(x$model, newdata = data, type = "response")
  }
  p <- p - x$shift
  return(p)
}

## pb2() using smooth.construct.
pb2 <- function(x, df = NULL, lambda = NULL, max.df = NULL, control = pb2.control(...), ...)
{
  m <- c(control$degree, control$order)
  k <- control$inter + m[1L] + 1L
  sx <- s(x, bs = "ps", m = m, k = k, sp = lambda)
  sx$control <- control
  sx$label <- gsub("s(", "pb2(", sx$label, fixed = TRUE)
  sx$orig.label <- sx$label
  sx$localML <- TRUE
  return(sx)
}

pb2.control <- function(inter = 20, degree = 3, order = 2,
  start = 10, method = c("ML", "GAIC", "GCV"), k = 2, ...) 
{
  if(inter <= 0) {
    warning("the value of inter supplied is less than 0, the value of 10 was used instead")
    inter <- 10
  }

  if(degree <= 0) {
    warning("the value of degree supplied is less than zero or negative the default value of 3 was used instead")
    degree <- 3
  }

  if(order < 0) {
    warning("the value of order supplied is zero or negative the default value of 2 was used instead")
    order <- 2
  }

  if(k <= 0) {
    warning("the value of GAIC/GCV penalty supplied is less than zero the default value of 2 was used instead")
    k <- 2
  }

  method <- match.arg(method)

  return(list(inter = inter, degree = degree, order = order, start = start, method = method, k = k))
}

