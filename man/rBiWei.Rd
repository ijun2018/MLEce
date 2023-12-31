% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLEce.R
\name{rBiWei}
\alias{rBiWei}
\title{Generating random data for the bivariate Weibull distribution.}
\usage{
rBiWei(n, paras)
}
\arguments{
\item{n}{number of observations.}

\item{paras}{parameters of bivariate Weibull distribution (alpha1, beta1, alpha2, beta2, delta).}
}
\value{
\code{rBiWei} generates random deviates.The length of generated data is determined by \code{"n"}
}
\description{
Generating random data for the bivariate Weibull distribution.
}
\details{
\code{rBiWei} generates random number data for bivariate Weibull distribution.
}
\examples{
datt = rBiWei(n=50, c(4,3,3,4,0.6))
}
