\name{foreca}
\alias{foreca}
\title{Performs Forecastable Component Analysis}
\usage{
foreca(series, method = "EM", ...)
}
\arguments{
  \item{series}{multivariate time series}

  \item{method}{what method should be used; currently only
  the EM-like algorithm from the original paper is
  available (see References): \code{method = "EM"} calls
  \code{\link{foreca.EM}}.}

  \item{...}{additional arguments passed to the individual
  methods (currently only \code{\link{foreca.EM}}).}
}
\value{
  A list with similar output as
  \code{\link[stats]{princomp}}. Signals are ordered from
  most to least forecastable.
}
\description{
  \code{foreca} performs Forecastable Component Analysis on
  a given multivariate and returns an object of class
  \code{foreca}.
}
\examples{
\dontrun{
XX <- diff(log(EuStockMarkets[-c(1:50), ])) * 100
plot(ts(XX))
ff <- foreca(XX[, 1:3], n.comp = 3)

summary(ff)
plot(ff)
}
}
