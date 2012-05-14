\name{ForeCA-package}
\alias{ForeCA-package}
\alias{ForeCA}
\docType{package}
\title{Forecastable Component Analysis}
\description{
ForeCA is a blind source separation (BSS) technique that finds an optimal linear
transformation of time series such that the transformed series is as forecastable as
possible.

The measure of forecastability \eqn{\Omega(x_t)} is based on the entropy of the spectral density of 
a time series: higher entropy means less forecastable, lower entropy is more forecastable. With this objective function it is possible to find optimal linear combinations \eqn{y_t = \mathbf{w}' \mathbf{X}_t}, such that \eqn{y_t} is easy to predict.
}
\details{
\tabular{ll}{
Package: \tab ForeCA\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2012-05-18\cr
License: \tab GPL\cr
}

In this current (first release) version \code{\link{Omega}} is the main function. It estimates
the forecastability for a univariate time series and returns a value between \eqn{0\%} and \eqn{100\%}: white noise has \eqn{0} forecastability, a perfect sinusoid is \eqn{100\%} predictable (however, empirically it reaches \eqn{100\%} only if the estimated spectrum has only one peak at some frequency \eqn{\omega_j} and \eqn{\widehat{f}(\omega_k) = 0} for all \eqn{k\neq j} ).

The BSS implementation of ForeCA will be added in a future version.
}
\author{
Author and maintainer: Georg M. Goerg <gmg@stat.cmu.edu>
}
\references{
Goerg, G.M. (2012). \dQuote{ForeCA: Forecastable Component Analysis}. In preparation for submission (available on http://www.arxiv.org soon).
}

\keyword{ package }

\examples{
# see examples for the function 'Omega'
}
