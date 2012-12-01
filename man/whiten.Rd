\name{whiten}
\alias{whiten}
\title{whitens multivariate data}
\usage{
whiten(data)
}
\arguments{
  \item{data}{\eqn{T \times n} array}
}
\value{
  A list with the whitened data, the transformation, and
  other useful quantities
}
\description{
  This function transforms a multivariate signal
  \eqn{\mathbf{X}} with mean \eqn{\mu_X} and covariance
  matrix \eqn{\Sigma_X} to a \emph{whitened} signal
  \eqn{\mathbf{U}} with mean \eqn{0} and \eqn{\Sigma_U =
  I_n}. That is it makes the signal contemporaneously
  uncorrelated.

  Same as projecting it on the principal components (PCs)
  and then scaling all PCs to unit variance.
}
\details{
  The whitening is achieved by pre-multiplying
  \eqn{\mathbf{X}} with the inverse of the \emph{square
  root} of \eqn{\Sigma_X} using the eigen-decomposition
  \deqn{ \Sigma_X = V' \Lambda V, } where \eqn{\Lambda} is
  an \eqn{n \times n} matrix with the eigenvalues
  \eqn{\lambda_1, \ldots, \lambda_n} in the diagonal.  The
  \emph{inverse square root} is defined as
  \eqn{\Sigma_X^{-1/2} = V' \Lambda^{-1/2}}, where
  \eqn{\Lambda^{-1/2} = diag(\lambda_1^{-1/2}, \ldots,
  \lambda_n^{-1/2})}.
}
\examples{
set.seed(1)
XX <- matrix(rnorm(100), ncol = 2)
cov(XX)
UU <- whiten(XX)$U
cov(UU)
}
\keyword{manip}
