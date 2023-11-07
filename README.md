
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MLEce

<!-- badges: start -->
<!-- badges: end -->

Welcome to MLEce, an R package for asymptotic efficient closed-form
estimators for multivariate distributions.

# Description

The goal of MLEce is to provide asymptotic efficient closed-form
estimators (MLEces) for three multivariate distributions (gamma, Weibull
and Dirichlet) whose maximum likelihood estimators (MLEs) are not in
closed forms. Closed-form estimators are strong consistent, and have the
similar asymptotic normal distribution like MLEs. But the calculation of
MLEces are much faster than the corresponding MLEs.

## Installation

In the R program, you can install the MLEce package like so:

``` r
install.packages(MLEce)
library(MLEce)
```

## Examples

There are basic examples which show you how to solve common problems.
For the bivariate gamma distribution, the following example is
presented.

``` r
library(MLEce)
#bivariate gamma distribution
 data_BiGam <-  rBiGam(100, c(1,4,5))
 res_BiGam  <- MLEce(data_BiGam, "BiGam")
 summary(res_BiGam)
```

For the bivariate Weibull distribution, a real data example is provided
as follows.

``` r
 #bivariate Weibull distribution
 data(airquality)
 air_data <- airquality[ ,3:4]
 air_data[ ,2] <- air_data[ ,2]*0.1
 est_BiWei <- MLEce(air_data, "BiWei")
 print(est_BiWei)
```

For the multivariate Dirichlet distribution, an example about the
comparison between the MLEce and MLE is given as follows.

``` r
 #multivariate Dirichlet distribution
 data_Diri <- LaplacesDemon::rdirichlet(80, c(3,4,1,3,4))
 benchMLEce(data_Diri, distname="Dirichlet", methods=c("MLEce","MLE"))
```

More commands of the MLEce with corresponding examples can be found in
the CRAN.

# Licensing

MLEce is licensed for public/open source use under the GNU General
Public License, Version 2 (GPL-2).
