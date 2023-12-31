% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLEce.R
\name{plot.MLEce}
\alias{plot.MLEce}
\title{Providing some plots for effective closed-form estimators}
\usage{
\method{plot}{MLEce}(
  x,
  which = c(1, 2, 3, 4),
  ask = prod(par("mfcol")) < length(which) && dev.interactive(),
  ...
)
}
\arguments{
\item{x}{an object of class "MLEce" made by the function \code{MLEce}.}

\item{which}{if a subset of the plots is required, specify a subset of 1:4.}

\item{ask}{logical; if TRUE, the user is asked before each plot.}

\item{...}{not used, but exists because of the compatibility.}
}
\description{
\code{plot} method for a class "MLEce".
}
\details{
The boxplot for given data is presented first with \code{which=1}. For \code{which=2}, a contour line is drawn by the probability density function of the estimated parameter based on effective closed-form estimators. In the counter plot, the x-axis is the first column of data and the y-axis is the second column of data. For \code{which=3}, a marginally fitted probability density plot is given for the first column of input data. And a fitted line is added for the efficient closed-form estimator. For \code{which=4}, is a marginally fitted probability density plot is given like the former one for the second column of input data. Note that, marginally fitted probability density plots in \code{which=3} and \code{which=4} present comparisons between efficient closed form estimators (MLEces) and correlation based method estimators (CMEs) for the bivariate Weibull distribution. Note that this \code{plot} commend is limited at the bivariate distributions.
}
\examples{
data(flood)
est_BiGam <- MLEce(flood, "BiGam")
plot(est_BiGam, c(3))
air_data <- airquality[ ,3:4]
air_data[ ,2] <- air_data[ ,2]*0.1
est_BiWei <- MLEce(air_data, "BiWei")
plot(est_BiWei)
data(fossil_pollen) 
fossil_data <- cbind(fossil_pollen[,1]/100,rowSums(fossil_pollen[,-1]/100))
est_fossil <- MLEce(fossil_data, "Dirichlet")
plot(est_fossil,c(2))
}
