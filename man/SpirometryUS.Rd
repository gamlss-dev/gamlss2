\name{SpirometryUS}

\alias{SpirometryUS}

\title{Spirometry Measurements from NHANES 2007--2012}

\description{
  Various spirometry measurements from the National Health and Nutrition
  Examination Survey (NHANES) 2007--2012 along with covariates providing
  demographics and basic body measurements.
}

\usage{data("SpirometryUS", package = "gamlss2")}

\format{
A data frame containing 16596 observations on 13 variables.
\describe{
  \item{fvc}{Numeric. Forced vital capacity (FVC) in liters, i.e., the
    volume of air that can forcibly be blown out after full inspiration.}
  \item{fev1}{Numeric. Forced expiratory volume in 1 second (FEV1) in
    liters, i.e., the volume of air that can forcibly be blown out in the
    first second, after full inspiration.}
  \item{ratio}{Numeric. Ratio of FEV1 to FVC.}
  \item{pef}{Numeric, peak expiratory flow (PEF) in liters per second,
    i.e., the maximal flow (or speed) achieved during the maximally forced
    expiration initiated at full inspiration.}
  \item{fef}{Numeric. Forced expiratory flow (FEF) in liters per second,
    i.e., the flow (or speed) of air coming out of the lung during the
    middle portion (25\% to 75\%) of a forced expiration.}
  \item{volume}{Numeric. Extrapolated volume.}
  \item{fet}{Numeric. Forced expiratory time (FET) in seconds,
    i.e., the length of the expiration.}
  \item{gender}{Factor. Binary gender information with levels
    \code{female} and \code{male}.}
  \item{age}{Numeric. Age in years (rounded to quarters).}
  \item{weight}{Numeric. Body weight in kilograms.}
  \item{height}{Numeric. Body height in centimeters.}
  \item{bmi}{Numeric. Body mass index in kilograms per meter-squared,
    rounded to 2 decimal places.}
  \item{ethnicity}{Factor. Self-reported race and ethnicity information
    with levels \code{white}, \code{black}, \code{mexican} American,
    other \code{hispanic}, and \code{other} (including multi-racial).}
  \item{surveyweight}{Numeric. Full sample two-year exam survey weight.}
}
}

\details{
  In order to establish lung function reference equations, Zavorsky (2025)
  studies the dependence of three spirometry measurements (FVC, FEV1, and
  the FEV1/FVC ratio) on age, adjusted for height and weight and separately
  for females and males. He intends to show that a simple normally-distributed
  model with (piecewise) linear mean equation and (piecewise) constant
  variance suffices for obtaining an adequate distributional fit from which
  the 5\% quantile can be obtained as the so-called lower limit of normal (LLN).
  Actually, his comparison with GAMLSS -- using flexible predictors for both
  mean and variance along with a Box-Cox-transformed normal distribution --
  shows that GAMLSS leads to a similar fit for the mean but a much better fit
  for the LLN.
  
  Zavorsky's (2025) analyses are based on a data set that he derived from the
  National Health and Nutrition Examination Survey (NHANES) in the United
  States 2007--2012. From the entire available data from
  \url{https://wwwn.cdc.gov/nchs/nhanes/} he included those observations
  which met or exceeded the technical acceptability of the measurements
  for forced expiratory volume in 1 second (FEV1) and forced vital capacity
  (FVC). The data are described in a short communication published in the
  \emph{Data in Brief} journal and the accompanying spreadsheet in CSV format
  (comma-separated values) is available from \emph{Mendeley Data}.
  
  The data comprises observations from NHANES' \dQuote{Examination Data},
  in particular in \dQuote{Spirometry -- Pre and Post-Bronchodilator} and
  \dQuote{Body Measures}, plus accompanying \dQuote{Demographics Data}.
  Respondent sequence numbers are stored as row names.
  See the variable descriptions above for more details. Basic information
  about spirometry can be found for example in the Wikipedia at
  \url{https://en.wikipedia.org/wiki/Spirometry}.
}

\source{
  Zavorsky GS (2024). \dQuote{Refined NHANES 2007--2012 Spirometry Dataset
    for the Comparison of Segmented (Piecewise) Linear Models to That of
    GAMLSS}, \emph{Mendeley Data}, V1. \doi{10.17632/dwjykg3xww.1}
}

\references{
  Zavorsky GS (2024). \dQuote{A Refined Spirometry Dataset for Comparing
    Segmented (Piecewise) Linear Models to that of GAMLSS}.
    \emph{Data in Brief}, \bold{57}, 111062.
    \doi{10.1016/j.dib.2024.111062}

  Zavorsky GS (2025). \dQuote{Debunking the GAMLSS Myth: Simplicity Reigns
    in Pulmonary Function Diagnostics}. \emph{Respiratory Medicine},
    \bold{236}, 107836. \doi{10.1016/j.rmed.2024.107836}
}

\examples{
data("SpirometryUS", package = "gamlss2")
summary(SpirometryUS)
}

\keyword{datasets}
