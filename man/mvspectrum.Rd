\name{mvspec2mvspectrum}
\alias{mvspec2mvspectrum}
\alias{mvspectrum}
\alias{normalize_mvspectrum}
\alias{SDF2mvspectrum}
\title{Estimates spectrum of multivariate time series and returns a 3D array}
\usage{
mvspec2mvspectrum(mvspec.output)

mvspectrum(series, method = "multitaper", normalize = FALSE, smoothing = FALSE, 
    Sigma_X = NULL, ...)

normalize_mvspectrum(mvspectrum.output, Sigma_X = NULL)

SDF2mvspectrum(sdf.output)
}
\arguments{
  \item{mvspec.output}{output from
  \code{\link[astsa]{mvspec}}}

  \item{series}{univariate or multivariate time series}

  \item{method}{method for spectrum estimation; see
  \code{method} argument in \code{\link[sapa]{SDF}} or use
  \code{'mvspec'} for estimation via
  \code{\link[astsa]{mvspec}} from the \pkg{astsa}
  package.}

  \item{normalize}{indicator; if \code{TRUE} the spectrum
  will be normalized such that it adds up to \eqn{0.5} (for
  univariate spectra) or to \eqn{0.5 \times I_n} (for
  multivariate spectra from \eqn{n} time series).}

  \item{smoothing}{indicator; if \code{TRUE} the spectrum
  gets additionally smoothed using a nonparametric smoother
  from \code{\link[mgcv]{gam}} with an automatically chosen
  (cross-validation) smoothing parameter.}

  \item{\dots}{additional arguments passed to
  \code{\link[sapa]{SDF}} or \code{\link[astsa]{mvspec}}
  (e.g. \code{taper})}

  \item{mvspectrum.output}{an object of class
  \code{"mvspectrum"}}

  \item{Sigma_X}{optional; covariance matrix of
  \code{series}}

  \item{sdf.output}{an object of class \code{"SDF"} from
  \code{\link[sapa]{SDF}}}
}
\value{
  An array of dimension \eqn{nfreqs \times nseries \times
  nseries}.

  \code{mvspectrum} returns a 3D array of dimension
  \eqn{nfreqs \times nseries \times nseries}, where
  \itemize{ \item nfreqs is the number of frequencies \item
  nseries are number of series (number of columns in
  \code{series}). }

  Returns a normalized spectrum. \describe{
  \item{univariate:}{output adds up to \eqn{0.5}.}
  \item{multivariate:}{output adds up to \eqn{0.5 \, I_n},
  where \eqn{I_n} is the \eqn{n \times n} identitiy
  matrix.} }
}
\description{
  \code{mvspec2mvspectrum} converts output from the
  \code{\link[astsa]{mvspec}} function to the
  \code{\link{mvspectrum}} output.

  Spectrum of a multivariate time series are matrix-valued
  functions of the frequency \eqn{\lambda}.
  \code{mvspectrum} estimates these spectra and puts them
  in a 3D array of dimension \eqn{nfreqs \times nseries
  \times nseries}.

  \code{SDF2mvspectrum} converts \code{\link[sapa]{SDF}}
  output to a 3D array with number of frequencies in the
  first dimension and the spectral density matrix in the
  other two (see Details below).
}
\details{
  The \code{\link[astsa]{mvspec}} function has frequencies
  in the first dimension, \code{\link{mvspectrum}} in the
  last.  \code{mvspec2mvspectrum} simply reshapes this
  \eqn{3D} array to an array of size \eqn{nfreqs \times
  nseries \times nseries}.

  The \code{\link[sapa]{SDF}} function only returns the
  upper diagonal (including diagonal), since spectrum
  matrices are Hermitian. For efficient vectorized
  computations, however, the full matrices are required.
  Thus \code{SDF2mvspectrum} converts SDF output to a 3D
  array with number of frequencies in the first dimension
  and the spectral density matrix in the other two (see
  Details below).

  \code{SDF2mvspectrum} is typically not called by the
  user, but by \code{\link{mvspectrum}}.
}
\examples{
# use SDF first and then SDF2mvspectrum
set.seed(1)
XX <- cbind(rnorm(100), arima.sim(n = 100, list(ar = 0.9)))
ss3d <- mvspectrum(XX)
dim(ss3d)

ss3d[2, , ]  # spectrum at omega_1; in general complex-valued, but Hermitian
identical(ss3d[2, , ], Conj(t(ss3d[2, , ])))  # is Hermitian
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
