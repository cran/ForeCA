\name{spectral.entropy}
\alias{spectral.entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimates the spectral entropy of a time series}
\description{
Estimates \emph{spectral entropy} from a univariate (or multivariate) time series. See Details below.
}
\usage{
spectral.entropy(series, spectrum_method = "wosa", base = exp(1), entropy_method = "MLE", threshold = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{ univariate time series (multivariate also supported but generally not used).}
  \item{spectrum_method}{ method for spectrum estimation; see \code{method} argument in \code{\link[sapa]{SDF}} 
  }
  \item{base}{base of the logarithm; default \code{base = exp(1)}.}
  \item{entropy_method}{ method to estimate the entropy from discrete probabilities \eqn{p_i}; here those 'probabilities' are the spectral density evaluated at the Fourier frequencies: \eqn{p_i = \widehat{f}(\omega_j)}.}
  \item{threshold}{ values of spectral density below \code{threshold} are set to \eqn{0}; default \code{threshold = 0}.}
  \item{\dots}{ optional additional arguments passed to \code{\link[sapa]{SDF}} }
}

\details{
The \emph{spectral entropy} is defined as the entropy of the spectral density \eqn{f_x(\lambda)} of a stationary process \eqn{x_t}: 
\deqn{
H_s(x_t) = - \int_{-\pi}^{\pi} f_x(\lambda) \log f_x(\lambda) d \lambda \geq 0.
}
}
\value{A non-negative real value for the spectral entropy \eqn{H_s(x_t)}. }

\references{ Jerry D. Gibson and Jaewoo Jung (2006). \dQuote{The Interpretation of Spectral Entropy Based Upon Rate Distortion Functions}. IEEE International Symposium on Information Theory, pp. 277-281. }

\author{
Georg M. Goerg
}

\seealso{
\code{Omega}
}
\examples{
set.seed(1)
eps = rnorm(100)
spectral.entropy(eps)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{univar}% __ONLY ONE__ keyword per line
