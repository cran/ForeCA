\name{foreca}
\alias{foreca}
\title{Runs Forecastable Component Analysis}
\usage{
foreca(series, algorithm.type = c("EM"), ...)
}
\arguments{
  \item{series}{multivariate time series of dimension
  \eqn{T \times K}}

  \item{algorithm.type}{specifies the algorithm that should
  be used to estimate forecastable components. Currently
  only the EM-like algorithm from the original Goerg (2013)
  paper is available (see References): \code{algorithm.type
  = "EM"} calls \code{\link{foreca.EM}}.}

  \item{...}{additional arguments passed to the available
  ForeCA algorithms.}
}
\value{
An object of class \code{foreca}: A list with similar
output as \code{\link[stats]{princomp}}. Signals are
ordered from most to least forecastable.
}
\description{
\code{foreca} performs Forecastable Component Analysis
(ForeCA) on a \eqn{K}-dimensional time series with \eqn{T}
observations.
}
\examples{
\dontrun{
XX <- diff(log(EuStockMarkets[-c(1:50),])) * 100
plot(ts(XX))
ff <- foreca(XX[,1:3], n.comp = 3)

summary(ff)
plot(ff)
}
}

