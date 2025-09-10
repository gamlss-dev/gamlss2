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
## K-fold cross-validation
cv_gamlss2(..., data, folds = 5,
  metric = log_pdf_metric, parallel = FALSE, simplify = TRUE)

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
  \item{data}{a \code{data.frame} containing the variables in the model.}
  \item{folds}{either an integer specifying the number of folds,
    or a list, matrix, or data frame of index sets for test folds.
    Defaults to 5.}
  \item{metric}{a function of the form \code{metric(model, data)} returning
    a score for the given fitted model and test data. Defaults to \code{log_pdf_metric}.}
  \item{parallel}{logical. If \code{TRUE}, computation is carried
    out in parallel using \pkg{future.apply}.}
  \item{simplify}{logical. If \code{TRUE}, results are returned in
    a simplified vector or data frame depending on the metric output.}
  \item{model}{a fitted \code{gamlss2} model.}
  \item{data}{a \code{data.frame} for evaluating predictions or residuals.}
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
\dontshow{ if(!requireNamespace("distributions3")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }}
\dontrun{## load the abdominal circumference data
data("abdom", package = "gamlss.data")

## cross-validation using the NO distribution
## only model the mean with s(x)
cv1 <- cv_gamlss2(y ~ s(x), data = abdom, family = NO)

## now, also model the standard deviation with s(x)
cv2 <- cv_gamlss2(y ~ s(x) | s(x), data = abdom, family = BCT)

## evaluate log-likelihood
sum(cv1$score)
sum(cv2$score)
}}

