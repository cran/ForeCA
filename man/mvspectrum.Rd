\name{SDF2mvspectrum}
\alias{SDF2mvspectrum}
\alias{mvspec2mvspectrum}
\alias{mvspectrum}
\alias{normalize_mvspectrum}
\title{Estimates spectrum of multivariate time series}
\usage{
SDF2mvspectrum(sdf.output)

mvspec2mvspectrum(mvspec.output)

mvspectrum(series, method = c("multitaper", "direct", "lag window", "wosa",
  "mvspec", "ar"), normalize = FALSE, smoothing = FALSE, Sigma.X = NULL,
  ...)

normalize_mvspectrum(mvspectrum.output, Sigma.X = NULL)
}
\arguments{
  \item{sdf.output}{an object of class \code{"SDF"} from
  \code{\link[sapa]{SDF}}}

  \item{mvspec.output}{output from
  \code{\link[astsa]{mvspec}}}

  \item{series}{univariate or multivariate time series}

  \item{method}{method for spectrum estimation; see
  \code{method} argument in \code{\link[sapa]{SDF}} or use
  \code{'mvspec'} for estimation via
  \code{\link[astsa]{mvspec}} from the \pkg{astsa}
  package.}

  \item{normalize}{logical; if \code{TRUE} the spectrum
  will be normalized such that it adds up to \eqn{0.5} (for
  univariate spectra) or to \eqn{0.5 \times I_n} (for
  multivariate spectra from \eqn{n} time series).}

  \item{smoothing}{logical; if \code{TRUE} the spectrum
  gets additionally smoothed using a nonparametric smoother
  from \code{\link[mgcv]{gam}} with an automatically chosen
  (cross-validation) smoothing parameter.}

  \item{\dots}{additional arguments passed to
  \code{\link[sapa]{SDF}} or \code{\link[astsa]{mvspec}}
  (e.g., \code{taper})}

  \item{mvspectrum.output}{an object of class
  \code{"mvspectrum"}}

  \item{Sigma.X}{optional; covariance matrix of
  \code{series}}
}
\value{
\code{mvspec2mvspectrum} returns the same object as
\code{\link{mvspectrum}}.

\code{mvspectrum} returns a 3D array of dimension
\eqn{num.freqs \times K \times K}, where \itemize{ \item
num.freqs is the number of frequencies \item K is the
number of series (columns in \code{series}). }

\code{normalize_mvspectrum} returns a normalized spectrum:
\describe{ \item{univariate:}{output adds up to \eqn{0.5}.}
\item{multivariate:}{output adds up to \eqn{0.5 \, I_n},
where \eqn{I_n} is the \eqn{n \times n} identity matrix.} }
}
\description{
\code{SDF2mvspectrum} converts \code{\link[sapa]{SDF}}
output to a 3D array with number of frequencies in the
first dimension and the spectral density matrix in the
other two (see Details below).

\code{mvspec2mvspectrum} converts output from the
\code{\link[astsa]{mvspec}} function to the
\code{\link{mvspectrum}} output.

Spectra of a multivariate time series are matrix-valued
functions of the frequency \eqn{\lambda}. \code{mvspectrum}
estimates these spectra and puts them in a 3D array of
dimension \eqn{num.freqs \times K \times K}.
}
\details{
The \code{\link[sapa]{SDF}} function only returns the upper
diagonal (including diagonal), since spectrum matrices are
Hermitian. For fast vectorized computations, however, the
full matrices are required. Thus \code{SDF2mvspectrum}
converts SDF output to a \eqn{3D} array with number of
frequencies in the first dimension and the spectral density
matrix in the latter two.

\code{SDF2mvspectrum} is typically not called by the user,
but by \code{\link{mvspectrum}}.

The \code{\link[astsa]{mvspec}} function returns the
multivariate spectrum in a \eqn{3D} array with frequencies
in the first dimension, whereas \code{\link{mvspectrum}} in
the last.  \code{mvspec2mvspectrum} simply reshapes the
former to the latter array.
}
\examples{
# use SDF first and then SDF2mvspectrum
set.seed(1)
XX <- cbind(rnorm(100), arima.sim(n = 100, list(ar= 0.9)))
ss3d <- mvspectrum(XX)
dim(ss3d)

ss3d[2,,] # at omega_1; in general complex-valued, but Hermitian
identical(ss3d[2,,], Conj(t(ss3d[2,,]))) # is Hermitian
xx <- rnorm(1000)
var(xx)
mean(mvspectrum(xx, normalize = FALSE, method = "direct")) * 2

mean(mvspectrum(xx, normalize = FALSE, method = "wosa")) * 2
mean(mvspectrum(xx, normalize = TRUE, method = "wosa")) * 2
}
\references{
See References in \code{\link[stats]{spectrum}},
\code{\link[sapa]{SDF}}, \code{\link[astsa]{mvspec}}.
}
\keyword{manip}
\keyword{ts}

