\name{re}
\alias{re}


\title{Random Effects}

\description{
There are two ways of fitting a random effect within \code{gamlss2}. The first, using  \code{s()}, is for a simple random effect, that is, when only one factor is entered the model as a smoother. This method uses the function \code{s()} of the package \pkg{mgcv} with argument \code{bs = "re"}. For example, if \code{area} is factor with several levels, \code{s(area, bs = "re")} will sringh the levels of \code{area} towards their mean level. The second, more general way, allows to fit more complicated random effect models using the function \code{re()}. The function \code{re()} is an interface connecting \code{gamlss2} with the specialised package for random effects \pkg{nlme}.

Here we document only the \code{re()} function only but we also give examples using \code{s(..., bs = "re")}. 
}

\usage{
re(fixed =~ 1, random = NULL, ...)
}

\arguments{
 
\item{fixed}{A formula that specifies the fixed effects of the \code{\link{nlme}{lme}} model. In most cases, this can also be included in the \code{\link{gamlss2}} parameter formula.}
\item{random}{A formula specifying the random effect part of the model,
    as in the \code{\link{nlme}{lme()}} function.}
\item{\dots}{ For the \code{re()} function, the dots argument is used
    to specify additional control arguments for the \code{\link{nlme}{lme}} function,
    such as the \code{method} and \code{correlation} arguments.}
}

\value{
  Function \code{\link[mgcv]{s}} with \code{bs = "re"} returns a smooth specification
  object of class \code{"re.smooth.spec"}, see
  also \code{\link[mgcv]{smooth.construct.re.smooth.spec}}.

  The \code{re()} function returns a special model term specification object,
  see \code{\link{specials}} for details.
}

\details{
  Both functions set up model terms that can be estimated using a backfitting algorithm, e.g.,
  the default \code{\link{RS}} algorithm.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link[mgcv]{smooth.construct.re.smooth.spec}}, \code{\link[mgcv]{s}}, \code{\link[nlme]{lme}}
 }

\examples{

\dontshow{ if(!requireNamespace("nlme")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}

## orthdontic measurement data
data("Orthodont", package = "nlme")

## model using lme()
m <- lme(distance ~ I(age-11), data = Orthodont,
  random =~ I(age-11) | Subject, method = "ML")

## using re(), function I() is not supported,
## please transform all variables in advance
Orthodont$age11  <- Orthodont$age - 11

## estimation using the re() constructor
b <- gamlss2(distance ~ s(age,k=3) + re(random =~ age11 | Subject),
  data = Orthodont)

## compare fitted values
plot(fitted(b, model = "mu"), fitted(m))
abline(0, 1, col = 4)

## extract summary for re() model term
st <- specials(b, model = "mu", elements = "model")
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
f <- follicles ~ s(Time) + re(random =~ sin1 | Mare,
  correlation = corARMA(q = 2), control = lmeControl(maxIter = 100))

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
}

\keyword{utilities}
\keyword{regression}

