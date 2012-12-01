\name{mvspectrum2wcov}
\alias{mvspectrum2wcov}
\title{Compute (weighted) covariance matrix from frequency spectrum}
\usage{
mvspectrum2wcov(mvspectrum.output, weights = NULL)
}
\arguments{
  \item{mvspectrum.output}{output of class
  \code{"mvspectrum"}}

  \item{weights}{a weighting sequence; if \code{NULL} equal
  weights will be used.}
}
\value{
  An \eqn{n \times n} matrix.
}
\description{
  Estimates the (weighted) covariance matrix from the
  frequency spectrum (see Details below).
}
\details{
  The covariance matrix of a multivariate time series
  equals the average of the power spectrum \deqn{ \Sigma_X
  = \int_{-\pi}^{\pi} S_{X}(\lambda) d \lambda. }

  We can therefore obtain a generalized covariance matrix
  by using a weighted average \deqn{ \tilde{\Sigma}_X =
  \int_{-\pi}^{\pi} K(\lambda) S_{X}(\lambda) d \lambda, }
  where \eqn{K(\lambda)} is a kernel symmetric around
  \eqn{0} which integrates to \eqn{1}. This allows us to
  filter out particular frequencies that we do (or don't)
  want to be part of the covariance matrix estimator.

  For ForeCA this function is especially important as we
  use \deqn{ K(\lambda) = -\log f_y(\lambda), } as the
  \emph{weights} (they don't add up to \eqn{1}!).
}
\examples{
# use SDF first and then SDF2mvspectrum
XX <- diff(log(EuStockMarkets)) * 100
set.seed(1)
XX <- cbind(rnorm(100), arima.sim(n = 100, list(ar = 0.9)))
SS <- mvspectrum(XX, "direct")

Sigma.hat <- cov(XX)
Sigma.hat.freq <- mvspectrum2wcov(SS)

image(Sigma.hat)
image(Sigma.hat.freq)
}
\seealso{
  \code{\link{mvspectrum}}
}
\keyword{ts}
