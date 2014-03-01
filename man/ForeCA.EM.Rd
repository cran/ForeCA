\name{foreca.EM}
\alias{foreca.EM}
\alias{foreca.EM.opt_weightvector}
\title{EM-like algorithm to estimate optimal ForeCA transformation}
\usage{
foreca.EM(series, spectrum.method = c("multitaper", "direct", "lag window",
  "wosa", "mvspec"), n.comp = 2, kernel = NULL, plot = TRUE,
  threshold = 0, control = list(tol = 1e-06, max.iter = 200, num.starts =
  10), smoothing = FALSE, ...)

foreca.EM.opt_weightvector(series, f.U = NULL,
  spectrum.method = c("multitaper", "direct", "lag window", "wosa", "mvspec"),
  entropy.method = "MLE", control = list(tol = 1e-06, max.iter = 100,
  num.starts = 10), kernel = NULL, threshold = 0, smoothing = FALSE, ...)
}
\arguments{
  \item{n.comp}{number of components to be extracted}

  \item{plot}{logical; if \code{TRUE} a plot of the current
  optimal solution \eqn{\mathbf{w}_i^*} will be shown and
  updated in each iteration of the EM-like algorithm.}

  \item{series}{a \eqn{T \times K} array containing a
  multivariate time series. Can be a \code{matrix},
  \code{data.frame}, or a multivariate \code{ts} object.}

  \item{f.U}{multivariate spectrum of class
  \code{mvspectrum} with \code{normalize = TRUE}}

  \item{spectrum.method}{string; method to estimate the
  spectrum; see \code{\link{mvspectrum}}.}

  \item{entropy.method}{string; method to estimate the
  entropy; see \code{\link{discrete_entropy}}.}

  \item{control}{a list with control parameters for the
  iterative EM algorithm: i) \code{tol} sets the tolerance
  level for convergence; ii) \code{max.iter} sets the
  maximum number of iterations; and iii) \code{num.starts}
  determines how many random starts should be done to avoid
  local minima (default:\code{5})}

  \item{kernel}{an R function \code{kernel =
  function(lambda) ...} that weighs the Fourier
  frequencies; if \code{NULL} (default) all frequencies get
  equal weight.}

  \item{threshold}{set estimate of spectral density below
  \code{threshold} to 0 (and renormalize so it adds up to
  one)}

  \item{smoothing}{logical; if \code{TRUE} the spectral
  density estimate will be smoothed using
  \code{\link[mgcv]{gam}}. See \code{\link{mvspectrum}} for
  further details.}

  \item{...}{other arguments passed to
  \code{\link{mvspectrum}}}
}
\value{
A list with similar output as
\code{\link[stats]{princomp}}. Signals are ordered from
most to least forecastable.
}
\description{
\code{foreca.EM} estimates the optimal transformation of
multivariate time \code{series} to obtain forecastable
signals using an EM-like algorithm.
}
\examples{
\dontrun{
XX <- diff(log(EuStockMarkets[-c(1:50),])) * 100
plot(ts(XX))
foreca <- foreca.EM(XX[,1:3], n.comp = 3)

summary(foreca)
plot(foreca)
}
\dontrun{
XX <- diff(log(EuStockMarkets)) * 100
one.weight <- foreca.EM.opt_weightvector(XX, smoothing = FALSE)

plot(one.weight)
}
}
\keyword{manip}

