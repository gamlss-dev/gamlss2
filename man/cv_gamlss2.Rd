\name{cv_gamlss2}
\alias{cv_gamlss2}
\alias{log_pdf_metric}
\alias{rqres_metric}
\alias{mse_metric}

\title{Cross Validation for \pkg{gamlss2} Models}

\description{
\code{cv_gamlss2()} implements K-fold cross validation for models fitted with
\code{\link{gamlss2}}. Different scoring rules can be supplied via the
\code{metric} argument. Convenience metric functions
(\code{log_pdf_metric()}, \code{rqres_metric()}, \code{mse_metric()}) are provided.
}

\usage{
## k-fold cross-validation
cv_gamlss2(..., data, folds = 5,
  metric = log_pdf_metric, simplify = TRUE, cores = 1L)

## log-pdf for each observation
log_pdf_metric(model, data)

## randomized quantile residuals
rqres_metric(model, data)

## mean squared error
mse_metric(model, data)
}

\arguments{
  \item{...}{model specification passed to \code{\link{gamlss2}} such as
  \code{formula}, \code{family}, etc.}
  \item{data}{a \code{data.frame} containing the variables in the model.
    For functions supplied to argument \code{metric}, a
    \code{data.frame} for evaluating predictions or residuals.}
  \item{folds}{either an integer specifying the number of folds,
    or a list, matrix, or data frame of index sets for test folds.
    Defaults to 5.}
  \item{metric}{a function of the form \code{metric(model, data)} returning
    a score for the given fitted model and test data. Defaults to \code{log_pdf_metric}.}
  \item{simplify}{logical. If \code{TRUE}, results are returned in
    a simplified vector or data frame depending on the metric output.}
  \item{cores}{Positive integer, number of workers for evaluating folds.
    The default \code{cores = 1L} runs sequentially; larger values use
    \code{\link[future.apply]{future_lapply}}
    with a temporary \code{\link[future]{multisession}} plan. The caller's full
    future plan is restored on exit, including on errors. Parallel evaluation
    requires the suggested packages \pkg{future} and \pkg{future.apply}.}
  \item{model}{a fitted \code{gamlss2} model.}
}

\details{
\code{cv_gamlss2()} splits the data into training and test folds.
For each fold the model is fitted on the training data, and the chosen 
\code{metric} is evaluated on the held-out test data. By default,
the scoring rule is the log predictive density (\code{log_pdf_metric}),
but other metrics can be used, such as randomized quantile residuals 
(\code{rqres_metric}) or mean squared error of the conditional mean 
(\code{mse_metric}).

The function returns either a list of fold-wise results or, if
\code{simplify = TRUE}, a named vector or data frame aligned with the original
observations.

Set the random seed before calling \code{cv_gamlss2()} to reproduce randomly
generated folds. Parallel evaluation also uses reproducible random-number
streams for stochastic fits or metrics.
}

\value{
If \code{simplify = TRUE} and the metric returns scalars,
a named numeric vector of fold scores is returned. Otherwise a
data frame with fold membership and scores per observation is returned.  
The convenience metrics return a numeric vector of scores or residuals.
}

\seealso{
\code{\link{gamlss2}}, \code{\link[distributions3]{log_pdf}}
}

\examples{
\dontrun{## air quality data
air <- subset(airquality, !is.na(Ozone))

## set the seed
set.seed(123)

## cross-validation using the GA distribution
## only model mu parameter with s(Temp)
cv1 <- cv_gamlss2(
  Ozone ~ s(Temp),
  data = air,
  family = GA, maxit = 100
)

## now, add effects
cv2 <- cv_gamlss2(
  Ozone ~ s(Temp) + s(Wind) | s(Temp),
  data = air,
  family = GA,
  folds = cv1$fold, maxit = 100
)

## evaluate log-likelihood
sum(cv1$score)
sum(cv2$score)

## evaluate folds with two workers, reusing the same fold assignments
cv_parallel <- cv_gamlss2(
  Ozone ~ s(Temp), data = air, family = GA,
  folds = cv1$fold, maxit = 100, cores = 2L
)
sum(cv_parallel$score)
}}

