% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Omega.R
\name{Omega}
\alias{Omega}
\title{Estimate forecastability of a time series}
\usage{
Omega(
  series = NULL,
  spectrum.control = list(),
  entropy.control = list(),
  mvspectrum.output = NULL
)
}
\arguments{
\item{series}{a univariate time series; if it is multivariate, then
\code{\link{Omega}} works component-wise (i.e., same as \code{apply(series, 2, Omega)}).}

\item{spectrum.control}{list; control settings for spectrum estimation. 
See \code{\link{complete_spectrum_control}} for details.}

\item{entropy.control}{list; control settings for entropy estimation.
See \code{\link{complete_entropy_control}} for details.}

\item{mvspectrum.output}{an object of class \code{"mvspectrum"} representing
the multivariate spectrum of \eqn{\mathbf{X}_t} (not necessarily \code{normalize}d).}
}
\value{
A real-value between \eqn{0} and \eqn{100} (\%). \eqn{0} means not
forecastable (white noise); \eqn{100} means perfectly forecastable (a
sinusoid).
}
\description{
An estimator for the forecastability \eqn{\Omega(x_t)} of a univariate time series \eqn{x_t}.
Currently it uses a discrete plug-in estimator given the empirical spectrum (periodogram).
}
\details{
The \emph{forecastability} of a stationary process \eqn{x_t} is defined as
(see References)

\deqn{
\Omega(x_t) = 1 - \frac{ - \int_{-\pi}^{\pi} f_x(\lambda) \log f_x(\lambda) d \lambda }{\log 2 \pi} \in [0, 1]
}
where \eqn{f_x(\lambda)} is the normalized spectral \emph{density} of \eqn{x_t}.
In particular \eqn{ \int_{-\pi}^{\pi} f_x(\lambda) d\lambda = 1}.


For white noise \eqn{\varepsilon_t} forecastability
\eqn{\Omega(\varepsilon_t) = 0}; for a sum of sinusoids it equals \eqn{100} \%.
However, empirically it reaches \eqn{100\%} only if the estimated spectrum has
exactly one peak at some \eqn{\omega_j} and \eqn{\widehat{f}(\omega_k) = 0}
 for all \eqn{k\neq j}.

In practice, a time series of length \code{T} has \eqn{T} Fourier frequencies
which represent a discrete
probability distribution.  Hence entropy of \eqn{f_x(\lambda)} must be
normalized by \eqn{\log T}, not by \eqn{\log 2 \pi}.

Also we can use several smoothing techniques to obtain a less variance estimate of
\eqn{f_x(\lambda)}.
}
\examples{

nn <- 100
eps <- rnorm(nn)  # white noise has Omega() = 0 in theory
Omega(eps, spectrum.control = list(method = "pgram"))
# smoothing makes it closer to 0
Omega(eps, spectrum.control = list(method = "mvspec"))

xx <- sin(seq_len(nn) * pi / 10)
Omega(xx, spectrum.control = list(method = "pgram"))
Omega(xx, entropy.control = list(threshold = 1/40))
Omega(xx, spectrum.control = list(method = "mvspec"),
      entropy.control = list(threshold = 1/20))

# an AR(1) with phi = 0.5
yy <- arima.sim(n = nn, model = list(ar = 0.5))
Omega(yy, spectrum.control = list(method = "mvspec"))

# an AR(1) with phi = 0.9 is more forecastable
yy <- arima.sim(n = nn, model = list(ar = 0.9))
Omega(yy, spectrum.control = list(method = "mvspec"))

}
\references{
Goerg, G. M. (2013). \dQuote{Forecastable Component
 Analysis}. Journal of Machine Learning Research (JMLR) W&CP 28 (2): 64-72, 2013.
 Available at \url{http://jmlr.org/proceedings/papers/v28/goerg13.html}.
}
\seealso{
\code{\link{spectral_entropy}}, \code{\link{discrete_entropy}},
\code{\link{continuous_entropy}}
}
\keyword{math}
\keyword{univar}
