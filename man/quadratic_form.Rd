\name{quadratic_form}
\alias{quadratic_form}
\title{Computes quadratic form x' A x}
\usage{
quadratic_form(mat, vec)
}
\arguments{
  \item{mat}{\eqn{n \times n} matrix (real or complex)}

  \item{vec}{\eqn{n \times 1} vector (real or complex)}
}
\value{
A real/complex value (\eqn{\mathbf{x}' \mathbf{A}
\mathbf{x}}).
}
\description{
Computes the quadratic form \eqn{\mathbf{x}' \mathbf{A}
\mathbf{x}} for an \eqn{n \times n} matrix \eqn{\mathbf{A}}
and an \eqn{n}-dimensional vector \eqn{\mathbf{x}}, i.e., a
wrapper for \code{t(x) \%*\% A \%*\% x}.

Works with real and complex valued matrices/vectors.
}
\examples{
set.seed(1)
AA <- matrix(1:4, ncol = 2)
bb <- matrix(rnorm(2))
t(bb) \%*\% AA \%*\% bb
quadratic_form(AA, bb)
}
\keyword{math}
\keyword{univar}

