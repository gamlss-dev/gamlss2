\name{re}
\alias{re}

\title{Random Effects}

\description{
  Random effects can be included in \code{gamlss2} models in two ways.  
  The first uses \code{s()} with \code{bs = "re"} for simple random effects,
  i.e., when a single factor is entered into the model as a smoother.
  This approach relies on the \code{s()} function from the \pkg{mgcv} package.
  For example, if \code{area} is a factor with several levels,
  then \code{s(area, bs = "re")} shrinks the level-specific effects of
  \code{area} towards their overall mean.  

  The second, more general approach uses the model term constructor \code{re()},
  which provides an interface to the specialised random-effects functionality
  in the \pkg{nlme} package. This allows fitting more complex random-effects structures.  

  This documentation focuses on the \code{re()} function, but we also provide examples
  using \code{s(..., bs = "re")}.
}

\usage{
re(random, correlation = NULL, ...)
}

\arguments{
  \item{random}{A formula specifying the random effect part of the model,
    as in the \code{\link[nlme]{lme}} function.}
  \item{correlation}{An optional correlation object, see \code{\link[nlme]{lme}}.}
  \item{\dots}{For the \code{re()} function, the dots argument is used
    to specify additional control arguments for the \code{\link[nlme]{lme}} function,
    such as the \code{method} and \code{correlation} arguments.}
}

\value{
  Function \code{\link[mgcv]{s}} with \code{bs = "re"} returns a smooth
  specification object of class \code{"re.smooth.spec"}, see
  also \code{\link[mgcv]{smooth.construct.re.smooth.spec}}.

  The \code{re()} function returns a special model term specification object,
  see \code{\link{specials}} for details.
}

\details{
  Both functions set up model terms that can be estimated using a
  backfitting algorithm, e.g., the default \code{\link{RS}} algorithm.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link[mgcv]{smooth.construct.re.smooth.spec}},
  \code{\link[mgcv]{s}}, \code{\link[nlme]{lme}}
}

\examples{\dontshow{ if(!requireNamespace("nlme")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
\dontrun{## orthdontic measurement data
data("Orthodont", package = "nlme")

## model using lme()
m <- lme(distance ~ I(age-11), data = Orthodont,
  random =~ I(age-11) | Subject, method = "ML")

## using re(), function I() is not supported,
## please transform all variables in advance
Orthodont$age11  <- Orthodont$age - 11

## estimation using the re() constructor,
## setup formula first
f <- distance ~ age11 + re(~age11|Subject)

## estimate model
b <- gamlss2(f, data = Orthodont)

## compare fitted values
plot(fitted(b, model = "mu"), fitted(m))
abline(0, 1, col = 4)

## extract summary for re() model term
st <- specials(b, model = "mu", term = "age11", elements = "model")
summary(st)

## random intercepts and slopes with s() using AIC
a <- gamlss2(distance ~ s(age,k=3) + s(Subject, bs = "re") + s(Subject, age11, bs = "re"),
  data = Orthodont)

## compare fitted values
plot(fitted(b, model = "mu"), fitted(m))
points(fitted(a, model = "mu"), fitted(m), col = 2)
abline(0, 1, col = 4)

## more complicated correlation structures.
data("Ovary", package = "nlme")

## ARMA model
m <- lme(follicles ~ sin(2 * pi * Time) + cos(2 * pi * Time), data = Ovary, 
  random = pdDiag(~sin(2*pi*Time)), correlation = corARMA(q = 2))

## now with gamlss2(), transform in advance
Ovary$sin1 <- sin(2 * pi * Ovary$Time)
Ovary$cos1 <- cos(2 * pi * Ovary$Time)

## model formula
f <- follicles ~ sin1 + cos1 +
  re(~ 1 | Mare) +
  re(~ sin1 - 1 | Mare, correlation = corARMA(q = 2))

## estimate model
b <- gamlss2(f, data = Ovary)

## smooth random effects
f <- follicles ~ ti(Time) + ti(Mare, bs = "re") + 
  ti(Mare, Time, bs = c("re", "cr"), k = c(11, 5))

g <- gamlss2(f, data = Ovary)

## compare fitted values
par(mfrow = n2mfrow(nlevels(Ovary$Mare)), mar = c(4, 4, 1, 1))

for(j in levels(Ovary$Mare)) {
  ds <- subset(Ovary, Mare == j)

  plot(follicles ~ Time, data = ds)

  f <- fitted(b, model = "mu")[Ovary$Mare == j]
  lines(f ~ ds$Time, col = 4, lwd = 2)

  f <- fitted(g, model = "mu")[Ovary$Mare == j]
  lines(f ~ ds$Time, col = 3, lwd = 2)

  f <- fitted(m)[Ovary$Mare == j]
  lines(f ~ ds$Time, col = 2, lwd = 2)
}

## simulated data
set.seed(1328)

n <- 10000
k <- 500

## generate random effects
f <- as.factor(sample(1:k, size = n, replace = TRUE))
ref <- rnorm(k, sd = 0.6)

## random effects only for sigma
y <- rBCT(n, mu = 10, sigma = exp(-1 + ref[f]), tau = 10)

## estimate model
b <- gamlss2(y ~ 1 | re(~ 1 | f), family = BCT)

## extract fitted random effects
cb <- coef(b, full = TRUE)
cb <- cb[grepl("re", names(cb), fixed = TRUE)]

## compare
plot(cb, ref, main = round(mean((cb - ref)^2), 2))
abline(0, 1, lwd = 2, col = 4)
}}

\keyword{utilities}
\keyword{regression}

