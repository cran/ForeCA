\name{mvspectrum2wcov}
\alias{mvspectrum2wcov}
\title{Compute (weighted) covariance matrix from frequency spectrum}
\usage{
mvspectrum2wcov(mvspectrum.output, weights = 1)
}
\arguments{
  \item{mvspectrum.output}{output of class
  \code{"mvspectrum"}}

  \item{weights}{a sequence of weights for each frequency.
  By default uses weights that average out to \code{1}.}
}
\value{
An \eqn{n \times n} matrix.
}
\description{
Estimates the (weighted) covariance matrix from the
frequency spectrum (see Details).
}
\details{
The covariance matrix of a multivariate time series equals
the average of the power spectrum \deqn{ \Sigma_{X} =
\int_{-\pi}^{\pi} S_{X}(\lambda) d \lambda. }

A generalized covariance matrix can be obtained using a
weighted average \deqn{ \tilde{\Sigma}_X =
\int_{-\pi}^{\pi} K(\lambda) S_{X}(\lambda) d \lambda, }
where \eqn{K(\lambda)} is a kernel symmetric around \eqn{0}
which integrates to \eqn{1}. This allows one to remove or
amplify specific frequencies in the covariance matrix
estimate.

For ForeCA \code{mvspectrum2wcov} is especially important
as we use \deqn{ K(\lambda) = -\log f_y(\lambda), } as the
\emph{weights} (their average is not \eqn{1}!).
}
\examples{
# use SDF first and then SDF2mvspectrum
set.seed(1)
nn <- 20
YY <- cbind(rnorm(nn), arima.sim(n = nn, list(ar = 0.9)), rnorm(nn))
XX <- YY \%*\% matrix(rnorm(9), ncol = 3)
XX <- scale(XX, scale = FALSE, center = TRUE)

# sample estimate of covariance matrix
Sigma.hat <- cov(XX)
dimnames(Sigma.hat) <- NULL

# using the frequency spectrum
SS <- mvspectrum(XX, "wosa")
Sigma.hat.freq <- mvspectrum2wcov(SS)

plot(c(Sigma.hat/Sigma.hat.freq))
abline(h = 1)

image(Sigma.hat)
image(Sigma.hat.freq)
image(Sigma.hat / Sigma.hat.freq)
}
\seealso{
\code{\link{mvspectrum}}
}
\keyword{ts}

