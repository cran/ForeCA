\name{spectral_entropy}
\alias{spectral_entropy}
\title{Estimates spectral entropy of a time series}
\usage{
spectral_entropy(series, spectrum.method = c("multitaper", "direct",
  "lag window", "wosa", "mvspec", "ar"), spectrum.estimate = NULL,
  base = NULL, entropy.method = c("MLE"), threshold = 0,
  smoothing = FALSE, ...)
}
\arguments{
  \item{series}{univariate time series of length \eqn{T}
  (multivariate also supported but generally not used).}

  \item{spectrum.method}{method for spectrum estimation;
  see \code{method} argument in \code{\link{mvspectrum}}}

  \item{base}{base of the logarithm in the entropy
  estimate; default \code{base = NULL}. In this case it is
  set internally to \code{base = T}, such that the entropy
  is bounded above by \eqn{1}.}

  \item{entropy.method}{method to estimate the entropy from
  discrete probabilities \eqn{p_i}; here 'probabilities'
  are the spectral density evaluated at the Fourier
  frequencies, \eqn{\widehat{p}_i =
  \widehat{f}(\omega_i)}.}

  \item{spectrum.estimate}{optional; one can directly
  provide an estimate of the spectrum}

  \item{threshold}{values of spectral density below
  \code{threshold} are set to \eqn{0}; default
  \code{threshold = 0}.}

  \item{smoothing}{logical; if \code{TRUE} the spectrum
  will be smoothed with a nonparametric estimate using
  \code{\link[mgcv]{gam}} with an automatically chosen
  (cross-validation) smoothing parameter.}

  \item{\dots}{optional additional arguments passed to
  \code{\link{mvspectrum}}}
}
\value{
A non-negative real value for the spectral entropy
\eqn{H_s(x_t)}.
}
\description{
Estimates \emph{spectral entropy} from a univariate (or
multivariate) time series.
}
\details{
The \emph{spectral entropy} equals the Shannon entropy of
the spectral density \eqn{f_x(\lambda)} of a stationary
process \eqn{x_t}: \deqn{ H_s(x_t) = - \int_{-\pi}^{\pi}
f_x(\lambda) \log f_x(\lambda) d \lambda \geq 0, } where
the density is normalized such that \eqn{\int_{-\pi}^{\pi}
f_x(\lambda) d \lambda = 1}.

An estimate can be obtained using the periodogram (from
\code{\link{mvspectrum}}); thus using discrete, and not
continuous entropy.
}
\examples{
set.seed(1)
eps <- rnorm(100)
spectral_entropy(eps)

phi.v <- seq(-0.95, 0.95, by = 0.1)
SE <- matrix(NA, ncol = 3, nrow = length(phi.v))


for (ii in seq_along(phi.v)) {
 xx.temp <- arima.sim(n = 1000, list(ar = phi.v[ii]))
 SE[ii, 1] <- spectral_entropy(xx.temp, spectrum.method = "direct")
 SE[ii, 2] <- spectral_entropy(xx.temp, spectrum.method = "multitaper")
 SE[ii, 3] <- spectral_entropy(xx.temp, spectrum.method = "wosa")
}

matplot(phi.v, SE, type = "l", col = 1:3)
legend("bottom", c("direct", "multitaper", "wosa"), lty = 1:3, col = 1:3)


# AR vs MA
SE.arma <- matrix(NA, ncol = 2, nrow = length(phi.v))
SE.arma[, 1] <- SE[, 2]

for (ii in seq_along(phi.v)){
 yy.temp <- arima.sim(n = 1000, list(ma = phi.v[ii]))
 SE.arma[ii, 2] <- spectral_entropy(yy.temp, spectrum.method = "multitaper")
}

matplot(phi.v, SE.arma, type = "l", col = 1:2, xlab = "parameter")
abline(v = 0)
legend("bottom", c("AR(1)", "MA(1)"), lty = 1:2, col = 1:2)
}
\references{
Jerry D. Gibson and Jaewoo Jung (2006). \dQuote{The
Interpretation of Spectral Entropy Based Upon Rate
Distortion Functions}. IEEE International Symposium on
Information Theory, pp. 277-281.

L. L. Campbell, \dQuote{Minimum coefficient rate for
stationary random processes}, Information and Control, vol.
3, no. 4, pp. 360 - 371, 1960.
}
\seealso{
\code{\link{Omega}}, \code{\link{discrete_entropy}}
}
\keyword{ts}
\keyword{univar}

