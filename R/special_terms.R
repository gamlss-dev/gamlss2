## Special conditional inference ctree constructor.
tree <- function(formula, ...)
{
  stopifnot(requireNamespace("rpart"))
  st <- list()
  ctr <- list(...)
  st$control <- do.call(rpart::rpart.control, ctr)
  st$formula <- formula
  st$term <- all.vars(formula)
  st$label <- paste0("tree(", paste0(gsub(" ", "", as.character(formula)), collapse = ""), ")")
  st$data <- model.frame(formula)
  class(st) <- c("special", "tree")
  return(st)
}

## ctree fitting function for the backfitting algorithm.
special_fit.tree <- function(x, z, w, y, eta, j, family, control, ...)
{
  f <- update(x$formula, response_z ~ .)
  x$data$response_z <- z
  x$data$w <- w
  rval <- list(
    "model" = rpart::rpart(formula = f, data = x$data, weights = w,
      control = x$control)
  )
  rval$fitted.values <- predict(rval$model)
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift
  frame <- rval$model$frame
  leaves <- frame$var == "<leaf>"
  size <- sum(leaves)
  edf <- 2 * size - 1
  rval$edf <- edf
  class(rval) <- "tree.fitted"
  return(rval)
}

## A ct predict method.
special_predict.tree.fitted <- function(x, data, ...)
{
  p <- predict(x$model, newdata = data, type = "vector")
  p <- p - x$shift
  return(p)
}

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

## ps2() & pb2() wrapper function using s().
ps <- pb <- function(x, k = 20, ...)
{
  sx <- s(x, bs = "ps", k = k, ...)
  sx$term <- deparse(substitute(x))
  sx$label <- paste0("pb2(", sx$term, ")")
  sx$control <- list("criterion" = "ml")
  sx$localML <- TRUE
  return(sx)
}

## The constructor function is used in the formula
## when calling gamlss2().
n <- function(formula, ...)
{
  stopifnot(requireNamespace("nnet"))

  ## List for setting up the special model term.
  st <- list()

  ## List of control arguments.
  ctr <- list(...)
  if(is.null(ctr$size))
    ctr$size <- 50
  if(is.null(ctr$maxit))
    ctr$maxit <- 1000
  if(is.null(ctr$decay))
    ctr$decay <- 0.1
  if(is.null(ctr$trace))
    ctr$trace <- FALSE
  if(is.null(ctr$MaxNWts))
    ctr$MaxNWts <- 10000
  if(is.null(ctr$scale))
    ctr$scale <- TRUE

  ## Put all information together.
  st$control <- ctr
  st$formula <- formula
  st$term <- all.vars(formula)
  st$label <- paste0("n(", paste0(gsub(" ", "", as.character(formula)), collapse = ""), ")")
  st$data <- model.frame(formula)

  ## Scale per default!
  if(ctr$scale) {
    sx <- list()
    for(j in colnames(st$data)) {
      if(!is.factor(st$data[[j]])) {
        sx[[j]] <- range(st$data[[j]])
        st$data[[j]] <- (st$data[[j]] - sx[[j]][1]) / diff(sx[[j]])
      }
    }
    st$scalex <- sx
  }

  ## Assign the "special" class and the new class "n".
  class(st) <- c("special", "n")

  return(st)
}

## Set up the special "n" model term fitting function
special_fit.n <- function(x, z, w, control, ...)
{
  ## Model formula needs to be updated.
  .fnns <- update(x$formula, response_z ~ .)

  ## Assign current working response.
  x$data$response_z <- z
  x$data$weights_w <- w

  ## Possible weights from last iteration.
  Wts <- list(...)$transfer$Wts

  ## Estimate model.
  nnc <- parse(text = paste0('nnet::nnet(formula = .fnns, data = x$data, weights = weights_w,',
      'size = x$control$size, maxit = x$control$maxit, decay = x$control$decay,',
      'trace = x$control$trace, MaxNWts = x$control$MaxNWts, linout = TRUE',
      if(!is.null(Wts)) ', Wts = Wts)' else ')'))

  rval <- list("model" = eval(nnc))

  ## Get the fitted.values.
  rval$fitted.values <- predict(rval$model)

  ## Transferring the weights for the next backfitting iteration.
  ## Note, "transfer" can be used to transfer anything from one
  ## iteration to the next.
  rval$transfer <- list("Wts" = rval$model$wts)

  ## Center fitted values.
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift

  ## Degrees of freedom.
  rval$edf <- length(coef(rval$model))

  ## Possible scaling.
  rval$scalex <- x$scalex

  ## Assign class for predict method.
  class(rval) <- "n.fitted"

  return(rval)
}

## Finally, the predict method.
special_predict.n.fitted <- function(x, data, se.fit = FALSE, ...)
{
  if(!is.null(x$scalex)) {
    for(j in names(x$scalex)) {
      data[[j]] <- (data[[j]] - x$scalex[[j]][1]) / diff(x$scalex[[j]])
    }
  }
  p <- predict(x$model, newdata = data, type = "raw")
  p <- p - x$shift
  if(se.fit)
    p <- data.frame("fit" = p)
  return(p)
}

## Linear model terms wrapper.
lin <- function(x, ..., ridge = FALSE)
{
  x <- lab <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
  f <- try(as.formula(x), silent = TRUE)
  is_f <- TRUE
  if(inherits(f, "try-error")) {
    v <- as.list(substitute(list(...)))[-1]
    if(length(v)) {
      lab <- paste0(x, ",", as.character(unlist(v)))
      v <- c(x, as.character(unlist(v)))
    } else {
      v <- x
    }
    f <- as.formula(paste("~", paste(v, collapse = "+")))
    is_f <- FALSE
  }
  if(ridge) {
    lab <- paste0(lab,",ridge=", ridge)
  }
  lab <- gsub(" ", "", lab)
  v <- all.vars(f)
  sx <- list("formula" = f, "term" = v,
    "label" = paste0("lin(", lab, ")"),
    "by" = "NA", "dim" = length(v))
  if(!ridge) {
    sx$sp <- 1e-10
  }
  class(sx) <- "lin.smooth.spec"
  return(sx)
}

smooth.construct.lin.smooth.spec <- function(object, data, knots)
{
  object$X <- model.matrix(object$formula, data = if (is.list(data)) 
    data[all.vars(reformulate(names(data))) %in% all.vars(object$formula)]
    else data)
  if(any(grepl("(Intercept)", colnames(object$X), fixed = TRUE))) {
    object$X <- object$X[, -1L, drop = FALSE]
  }
  object$bs.dim <- ncol(object$X)
  object$S <- list(diag(object$bs.dim))
  object$rank <- object$bs.dim
  object$null.space.dim <- 0
  object$C <- matrix(0, 0, ncol(object$X))
  object$side.constrain <- FALSE
  class(object) <- "lin.effect"
  return(object)
}

Predict.matrix.lin.effect <- function(object, data) 
{
  if(is.list(data)) 
    data <- data[all.vars(reformulate(names(data))) %in% all.vars(object$formula)]
  X <- model.matrix(object$formula, model.frame(object$formula, data, na.action = na.pass))
  if(any(grepl("(Intercept)", colnames(X), fixed = TRUE))) {
    X <- X[, -1L, drop = FALSE]
  }
  X[!is.finite(X)] <- 0
  return(X)
}

## Random effects.
re <- function(fixed = ~ 1, random = NULL, ...)
{
  stopifnot(requireNamespace("nlme"))

  call <- match.call(expand.dots = FALSE)

  ## List for setting up the special model term.
  st <- list()

  ## List of control arguments.
  ctr <- list(...)
  if(is.null(ctr$method))
    ctr$method <- "ML"

  ## Put all information together.
  st$control <- ctr$control
  st$fixed <- fixed
  st$random <- random
  st$term <- c(all.vars(fixed), all.vars(random))
  st$label <- gsub(" ", "", deparse(call))
  st$label <- gsub("re(", "re2(", st$label, fixed = TRUE)
  st$method <- ctr$method
  ctr$method <- NULL
  st$correlation <- ctr$correlation
  ctr$correlation <- NULL

  ## Assign data.
  rm_v <- function(formula) {
    env <- environment(formula)
    formula <- formula(as.Formula(formula), drop = TRUE, collapse = TRUE)
    environment(formula) <- env
    return(formula)
  }

  df <- model.frame(rm_v(fixed))
  if(inherits(random, "formula"))
    dr <- model.frame(rm_v(random))
  else
    dr <- data.frame()

  if(nrow(df) < 1 && nrow(dr) > 0)
    st$data <- dr
  if(nrow(df) > 0 && nrow(dr) < 1)
    st$data <- df
  if(nrow(df) > 0 && nrow(dr) > 0)
    st$data <- cbind(df, dr)

  ## Assign the "special" class and the new class "n".
  class(st) <- c("special", "re")

  return(st)
}

## Set up the special "re" model term fitting function
special_fit.re <- function(x, z, w, control, ...)
{
  ## Model formula needs to be updated.
  .fnns <- update(x$fixed, response_z ~ .)

  ## Assign current working response.
  x$data$response_z <- z
  x$data$weights_w <- w

  ## Estimate model.
  rem <- parse(text = paste0(
    'nlme::lme(fixed = .fnns, data = x$data, random = x$random, weights = varFixed(~weights_w),',
    'method="', x$method, '",control=x$control,correlation=x$correlation,keep.data=FALSE)'))

  rval <- list("model" = eval(rem))

  ## Get the fitted.values.
  rval$fitted.values <- fitted(rval$model)

  ## Degrees of freedom.
  N <- sum(w != 0)
  rval$edf <- sum(w * (z - rval$fitted.values)^2)/(rval$model$sigma^2)

  ## Assign class for predict method.
  class(rval) <- "re.fitted"

  return(rval)
}

## The re() predict method.
special_predict.re.fitted <- function(x, data, se.fit = FALSE, ...)
{
  p <- predict(x$model, newdata = data, level = 0)
  if(se.fit)
    p <- data.frame("fit" = p)
  return(p)
}

