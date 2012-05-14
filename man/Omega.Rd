\name{Omega}
\alias{Omega}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plug-in estimator for measure of forecastability of a time series}
\description{
A plug-in estimator for \eqn{\Omega(x_t)} for a univariate time series \eqn{x_t}.

See References for details.
}
\usage{
Omega(series, spectrum_method = "multitaper", entropy_method = "MLE", threshold = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{ a univariate time series; if multivariate \code{Omega} works component-wise (i.e. same as using \code{apply(series, 2, Omega)}).}
  \item{spectrum_method}{ method to estimate the spectrum; see \code{method} argument in \code{\link[sapa]{SDF}}.}
  \item{entropy_method}{ method to estimate the entropy; see \code{method} argument in \code{discrete.entropy} for details.}
  \item{threshold}{ threshold for the entropy estimator; see \code{discrete.entropy} for details. }
  \item{\dots}{ optional additional arguments passed to \code{spectral.entropy} }
}

\value{
A real-value between \eqn{0} and \eqn{100} (\%). \eqn{0} means not forecastable (white noise); \eqn{100} means perfectly forecastable (a sinusoid).
}
\references{
Goerg, G.M. (2012). \dQuote{ForeCA: Forecastable Component Analysis}. In preparation for submission (available on http://www.arxiv.org soon).
}
\author{
Georg M. Goerg
}

\seealso{
\code{spectral.entropy}, \code{discrete.entropy}
}
\examples{
set.seed(1)
nn = 100
eps = rnorm(nn)
Omega(eps)
Omega(eps, spectrum_method = "wosa") # smoothing makes it closer to 0

xx = sin(1:nn * pi / 10)
Omega(xx) # direct (no smoothing of spectrum)
Omega(xx, threshold = 1/4) # direct (no smoothing of spectrum); and thresholding to single out the one important frequency
Omega(xx, spectrum_method = "wosa") # smoothing
Omega(xx, spectrum_method = "multitaper") # multitaper smoothing

# an AR(1) with phi = 0.9
yy = arima.sim(n = nn, model = list(ar=0.9))
Omega(yy, spectrum_method = "wosa")
Omega(yy, spectrum_method = "multitaper")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
\keyword{univar}% __ONLY ONE__ keyword per line
