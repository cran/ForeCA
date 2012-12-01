\name{ForeCA.EM}
\alias{ForeCA.EM}
\alias{ForeCA.EM.one_weightvector}
\title{EM-like Algorithm to Estimate Optimal Transformations via ForeCA}
\usage{
ForeCA.EM(series, spectrum_method = "multitaper", n.comp = 2, tol = 1e-04, max_iter = 200, 
    kernel = NULL, nstart = 10, plot = TRUE, threshold = 0, smoothing = FALSE, 
    ...)

ForeCA.EM.one_weightvector(series, f_U = NULL, spectrum_method = "multitaper", 
    entropy_method = "MLE", tol = 1e-04, max_iter = 100, kernel = NULL, nstart = 5, 
    threshold = 0, smoothing = FALSE, ...)
}
\arguments{
  \item{n.comp}{number of components to be extracted}

  \item{plot}{indicator; if \code{TRUE} a plot of the
  current optimal solution \eqn{\mathbf{w}_i^*} will be
  shown and the plot is updated for each new optimal
  weightvector.}

  \item{series}{an \eqn{T \times n} array containg a
  multivariate time series. Can be a \code{matrix},
  \code{data.frame}, or a multivariate \code{ts} object}

  \item{f_U}{multivariate spectrum of class
  \code{mvspectrum} with \code{normalize = TRUE}}

  \item{spectrum_method}{string; method to estimate the
  spectrum. See \code{\link{mvspectrum}}.}

  \item{entropy_method}{string; method to estimate the
  entropy. See \code{\link{discrete_entropy}}.}

  \item{tol}{tolerance level for convergence}

  \item{max_iter}{maximum number of iterations}

  \item{kernel}{an R function \code{kernel =
  function(lambda) ...} that weights the frequencies; if
  \code{NULL} (default) all frequencies get equal weight.}

  \item{nstart}{how many random starts should be done to
  avoid local minima?}

  \item{threshold}{set spectral density values below
  \code{threshold} to 0}

  \item{smoothing}{indicator; if \code{TRUE} then spectrum
  will be additionally smoothed using
  \code{\link[mgcv]{gam}}. See \code{\link{mvspectrum}} for
  details.}

  \item{...}{other arguments passed to
  \code{\link{mvspectrum}}}
}
\value{
  A list with similar output as
  \code{\link[stats]{princomp}}. Signals are ordered from
  most to least forecastable.
}
\description{
  \code{ForeCA.EM} estimates the optimal transformations to
  obtain forecastable signals from multivariate time series
  \code{series}
}
\examples{
\dontrun{
XX <- diff(log(EuStockMarkets[-c(1:50), ])) * 100
plot(ts(XX))
foreca <- ForeCA.EM(XX[, 1:3], n.comp = 3)

summary(foreca)
plot(foreca)
}
\dontrun{
XX <- diff(log(EuStockMarkets)) * 100
one_weight <- ForeCA.EM.one_weightvector(XX, smoothing = FALSE)

plot(one_weight)
}
}
\keyword{manip}
