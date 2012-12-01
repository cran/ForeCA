#'@title Estimate forecastability of a time series
#'
#'@description
#'A plug-in estimator for the forecastability \eqn{\Omega(x_t)} of 
#'a univariate time series \eqn{x_t}.
#'
#'@details
#'The \emph{forecastability} of a stationary process \eqn{x_t} is defined as 
#'(see References)
#'
#'\deqn{ 
#'\Omega(x_t) = 1 - \frac{ - \int_{-\pi}^{\pi} f_x(\lambda) \log f_x(\lambda) d \lambda }{\log 2 \pi} \in [0, 1]
#'}
#'where \eqn{f_x(\lambda)} is the spectral density of \eqn{x_t}. In particular
#'\eqn{ \int_{-\pi}^{\pi} f_x(\lambda) d\lambda = 1}.
#' 
#'For white noise \eqn{\varepsilon_t \sim WN(0, \sigma^2)} forecastability 
#'\eqn{\Omega(\varepsilon_t) = 0}; for a sum of sinusoids it equals \eqn{100} \%.
#'
#'However, empirically it reaches \eqn{100\%} only if the estimated spectrum has
#'exactly one peak at some \eqn{\omega_j} and \eqn{\widehat{f}(\omega_k) = 0}
#' for all \eqn{k\neq j}.
#'
#'In practice, we have \eqn{T} Fourier frequencies which represent a discrete 
#'probability distribution.  Hence entropy of \eqn{f_x(\lambda)} must be 
#'normalized by \eqn{\log T}, not by \eqn{\log 2 \pi}.
#'
#'@param series a univariate time series; if multivariate \code{\link{Omega}} works
#'component-wise (i.e. same as using \code{apply(series, 2, Omega)}).
#'@param spectrum_method method to estimate the spectrum; see \code{method}
#'argument in \code{\link{mvspectrum}}.
#'@param entropy_method method to estimate the entropy; see \code{method}
#'argument in \code{\link{discrete_entropy}} for details.
#'@param spectrum_estimate optionally one can directly pass an estimate of the
#'spectrum (rather than estimating it from \code{series} within
#' \code{\link{Omega}}).
#'@param threshold threshold for the entropy estimator; see
#'\code{\link{discrete_entropy}} for details.
#'@param smoothing indicator - see \code{\link{mvspectrum}}.
#'@param \dots optional arguments passed to \code{\link{spectral_entropy}}
#'@return A real-value between \eqn{0} and \eqn{100} (\%). \eqn{0} means not
#'forecastable (white noise); \eqn{100} means perfectly forecastable (a
#'sinusoid).
#'@export
#'@seealso \code{\link{spectral_entropy}}, \code{\link{discrete_entropy}}, 
#'\code{\link{continuous_entropy}}
#'@references Goerg, G.M. (2012). \dQuote{ForeCA: Forecastable Component
#'Analysis}. In preparation for submission (\url{arxiv.org/abs/1205.4591}).
#'@keywords math univar
#'@examples
#'
#'set.seed(1)
#'nn = 100
#'eps = rnorm(nn)
#'Omega(eps) # default is direct estimation; no smoothing
#'Omega(eps, spectrum_method = "wosa") # smoothing makes it closer to 0
#'
#'xx = sin(1:nn * pi / 10)
#'Omega(xx) # direct (no smoothing of spectrum)
#'Omega(xx, threshold = 1/40) # direct (no smoothing of spectrum); and thresholding to single out the one important frequency

#'Omega(xx, spectrum_method = "wosa") # smoothing
#'Omega(xx, spectrum_method = "wosa", threshold = 1/20) # smoothing and treshold
#'Omega(xx, spectrum_method = "multitaper") # multitaper smoothing
#'
#'set.seed(1)
#'# an AR(1) with phi = 0.5
#'yy = arima.sim(n = nn, model = list(ar=0.5))
#'Omega(yy, spectrum_method = "wosa")
#'Omega(yy, spectrum_method = "multitaper")
#'
#'# an AR(1) with phi = 0.9 is more forecastable
#'yy = arima.sim(n = nn, model = list(ar=0.9))
#'Omega(yy, spectrum_method = "wosa")
#'Omega(yy, spectrum_method = "multitaper")
#'
Omega <- function(series, 
                  spectrum_method = "direct", 
                  entropy_method = "MLE",
                  spectrum_estimate = NULL,
                  threshold = 0,
                  smoothing = FALSE, 
                  ...) {
  
  if (is.null(spectrum_estimate)) {
    series <- as.matrix(series)
    # scale to zero-mean and unit-variance
    series <- scale(series)
    nseries <- ncol(series)
    TT <- nrow(series)
  } else {
    TT <- nrow(as.matrix(spectrum_estimate))
    nseries <- ncol(as.matrix(spectrum_estimate))
  }
  
  if (is.null(threshold)) {
    threshold <- TT^(-1.5)
  }
  
  if (nseries > 1) {
    OMEGAs <- apply(series, 2, Omega, spectrum_method = spectrum_method, entropy_method = entropy_method, spectrum_estimate = spectrum_estimate, threshold = threshold, smoothing = smoothing)
    attr(OMEGAs, "unit") <- "%"
    return(OMEGAs)
  } else {
    H_s <- spectral_entropy(series, spectrum_method = spectrum_method, spectrum_estimate = spectrum_estimate, 
                            base = 2, entropy_method = entropy_method, threshold = threshold, smoothing = smoothing, 
                            ...)
    if (entropy_method == "MLE") {
      OMEGA <- 1 - H_s/(log(attr(H_s, "nfrequencies"), base = 2))
    }
    OMEGA <- c(OMEGA) * 100
    attr(OMEGA, "unit") <- "%"
    return(OMEGA)
  }
}
