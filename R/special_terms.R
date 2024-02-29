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

## pb2() using smooth.construct.
pb2 <- function(x, df = NULL, lambda = NULL, max.df = NULL, control = pb2.control(...), ...)
{
  m <- c(control$degree, control$order)
  k <- control$inter + m[1L] + 1L
  sx <- s(x, bs = "ps", m = m, k = k, sp = lambda)
  sx$term <- deparse(substitute(x))
  sx$control <- control
  sx$label <- paste0("pb2(", sx$term, ")")
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
    ctr$size <- 100
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
    for(j in colnames(data)) {
      if(!is.factor(data[[j]])) {
        sx[[j]] <- range(data[[j]])
        data[[j]] <- (data[[j]] - sx[[j]][1]) / diff(sx[[j]])
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
  rval$edf <- x$control$size

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

