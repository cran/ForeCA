\docType{package}
\name{ForeCA-package}
\alias{ForeCA}
\alias{ForeCA-package}
\title{Implementation of Forecastable Component Analysis (ForeCA)}
\description{
  ForeCA is a novel dimension reduction (DR) technique for
  multivariate time series.  ForeCA finds a linar
  combination \eqn{y_t = \mathbf{w}' \mathbf{X}_t} that is
  easy to forecast.  The measure of forecastability
  \eqn{\Omega(x_t)} (\code{\link{Omega}}) is based on the
  entropy of the spectral density \eqn{f_y(\lambda)} of
  \eqn{y_t}: higher entropy means less forecastable, lower
  entropy is more forecastable.

  The main function \code{\link{foreca}} runs ForeCA on a
  multivariate time series to find the most forecastable
  signals.

  Note that even though this version of the package has
  most functionality, some function naming conventions
  might change in future version.  In particular,
  underscores and ``.'' in the ``foreca.EM.x_y_z'' function
  names might change in future versions.  Please consult
  the manual of each new version.
}
\examples{
XX <- ts(diff(log(EuStockMarkets))[-c(1:1000), ])
Omega(XX)

plot(log(lynx, 10))
Omega(log(lynx, 10), spectrum_method = "direct")

\dontrun{
ff <- foreca(XX, n.comp = 2)
plot(ff)
summary(ff)
}
}
\author{
  Author and maintainer: Georg M. Goerg <gmg@stat.cmu.edu>
}
\references{
  Goerg, G.M. (2012). \dQuote{Forecastable Component
  Analysis}. Submitted for publication
  (\url{arxiv.org/abs/1205.4591}).
}
\keyword{package}
