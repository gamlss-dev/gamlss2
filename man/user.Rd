\name{user}
\alias{user}

\title{Reserved User Special Model Term}

\description{
  The model term constructor \code{user()} in \code{\link{gamlss2}} formulae
  may be used to set up user-defined model terms. See the example below.
}

\value{
  The value is determined by the user-defined constructor and
  should inherit from \code{c("special", "user")}.
}

\seealso{
  \code{\link{special_terms}}
}

\examples{\dontrun{
## a custom weighted local smoother
user <- function(formula, span = 0.75, degree = 2) {
  if(!inherits(formula, "formula"))
    formula <- as.formula(paste("~", deparse1(substitute(formula))))
  st <- list(
    formula = formula,
    data = model.frame(formula),
    control = list(span = span, degree = degree),
    term = all.vars(formula),
    label = paste0("user(", paste(all.vars(formula), collapse = "+"), ")")
  )
  class(st) <- c("special", "user")
  st
}

special_fit.user <- function(x, z, w, control, ...) {
  fit_formula <- update(x$formula, response_z ~ .)
  dat <- x$data
  dat$response_z <- z
  dat$weights_w <- w
  fit <- loess(
    fit_formula,
    data = dat,
    weights = weights_w,
    span = x$control$span,
    degree = x$control$degree
  )
  fitted <- as.numeric(predict(fit, newdata = dat))
  shift <- mean(fitted)
  out <- list(
    model = fit,
    fitted.values = fitted - shift,
    shift = shift,
    edf = fit$trace.hat
  )
  class(out) <- "user.fitted"
  out
}

special_predict.user.fitted <- function(x, data, se.fit = FALSE, ...) {
  fit <- as.numeric(predict(x$model, newdata = data)) - x$shift
  if(se.fit) return(data.frame(fit = fit))
  fit
}

air <- subset(airquality, !is.na(Ozone))

m <- gamlss2(
  Ozone ~ user(~ Temp, span = 0.7) | s(Temp),
  data = air,
  family = GA,
  trace = FALSE
)

plot(m)
}}

\keyword{models}

