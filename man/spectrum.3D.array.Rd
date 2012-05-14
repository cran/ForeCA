\name{spectrum.3D.array}
\alias{spectrum.3D.array}
\alias{SDF2spectrum.3D.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimates spectrum of multivariate time series and returns a 3D array}
\description{
Multivariate spectra are matrix-valued functions of the frequency \eqn{\lambda}. The \code{\link[sapa]{SDF}} function only returns the upper diagonal (including diagonal), since the matrices are Hermitian. For efficient vectorized computations, however, the full matrices are necessary; thus \code{SDF2spectrum.3D.array} converts SDF output to a 3D array with number of frequencies in the first dimension and the spectral density matrix in the other two (see Details below).

\code{SDF2spectrum.3D.array} is typically not called by the user, but by \code{\link{spectrum.3D.array}} which users should call instead directly: \code{spectrum.3D.array} estimates the spectrum of a (multivariate) time series and puts it in a 3D array. Basically a wrapper around \code{\link[sapa]{SDF}} and \code{SDF2spectrum.3D.array}.
}
\usage{
spectrum.3D.array(series, method = "wosa", frequency_0 = TRUE, ...)
SDF2spectrum.3D.array(sdf.output, frequency_0 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{univariate or multivariate time series }
  \item{method}{method for spectrum estimation; see \code{method} argument in \code{\link[sapa]{SDF}} }
  \item{frequency_0}{ should frequency \eqn{\omega_0 = 0} be included; default \code{TRUE}.}
  \item{sdf.output}{ output from \code{\link[sapa]{SDF}} }
  \item{\dots}{additional arguments passed to \code{\link[sapa]{SDF}}}
}

\value{
A 3D array of dimension \eqn{T \times nseries \times nseries}, where
  \item{T}{are the number of frequencies (if \code{frequency_0 = TRUE} then it is \eqn{T+1}).}
  \item{nseries}{are number of series (number of columns in \code{series}.}
%A 3D array containing the multivariate spectrum with number of frequencies
%in the first dimension and the spectral density matrix in the other two.
}
\references{
See References in \code{\link[stats]{spectrum}}.
}
\author{
Georg M. Goerg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[sapa]{SDF}}
}
\examples{
# use SDF first and then SDF2spectrum.3D.array
set.seed(1)
XX = cbind(rnorm(100), arima.sim(n = 100, list(ar= 0.9)))
ss = SDF(XX)
ss3d = SDF2spectrum.3D.array(ss)

dim(ss3d)

ss3d[1,,] # spectrum at first frequency (omega_0); real-valued for omega_0
ss3d[2,,] # spectrum at omega_1; in general complex-valued, but Hermitian

# or you can use spectrum.3D.array directly
ss3d_direct = spectrum.3D.array(XX)
dim(ss3d_direct)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{multivar}% __ONLY ONE__ keyword per line
