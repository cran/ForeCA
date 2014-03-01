\name{continuous_entropy}
\alias{continuous_entropy}
\title{Shannon entropy for a continuous pdf}
\usage{
continuous_entropy(pdf, a = NULL, b = NULL, base = 2)
}
\arguments{
  \item{pdf}{R function (for example, \code{function(x)
  return(dunif(x, 0, 1))}) for the pdf \eqn{p(x)} of a RV
  \eqn{X \sim p(x)}. Note that this function must integrate
  to \eqn{1} over the interval [\code{a}, \code{b}].}

  \item{base}{logarithm base; entropy is measured in
  ``nats'' for \code{base = e}; in ``bits'' if \code{base =
  2} (default).}

  \item{a}{lower integration limit}

  \item{b}{upper integration limit}
}
\value{
Entropy (real).
}
\description{
Computes the Shannon entropy \eqn{\mathcal{H}(p)} for a
continuous probability density function (pdf) \eqn{p(x)}.
}
\details{
The Shannon entropy of a continuous RV \eqn{X \sim p(x)} is
defined as \deqn{ \mathcal{H}(p) = -\int_{-\infty}^{\infty}
p(x) \log p(x) d x }

Contrary to the entropy of a discrete RV, continuous RVs
can have negative entropy (see Examples).
}
\examples{
my_density = function(x){
 dunif(x, -1, 1)
}

\dontrun{
# error since limits of integration were not specified
continuous_entropy(my_density)

# error since density does not integrate to '1' over [-0.5, 1]
continuous_entropy(my_density, -0.5, 1)
}
continuous_entropy(my_density, -1, 1) # = log(b - a)

# entropy of U(a,b) = log(b-a). Thus not necessarily positive anymore, e.g.
continuous_entropy(function(x) dunif(x, 0, 0.5), 0, 0.5) # log2(0.5)
}
\seealso{
\code{\link{discrete_entropy}}
}
\keyword{math}
\keyword{univar}

