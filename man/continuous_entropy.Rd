\name{continuous_entropy}
\alias{continuous_entropy}
\title{Shannon entropy for a continuous pdf}
\usage{
continuous_entropy(pdf = function(x) return(dunif(x, 0, 1)), a = NULL, b = NULL, base = 2)
}
\arguments{
  \item{pdf}{R function (\code{function(x) ...}) for the
  pdf \eqn{p(x)} of a RV \eqn{X \sim p(x)}}

  \item{base}{logarithm base; default \code{base = 2}.
  Entropy is measured in ``nats'' for \code{base = e}; in
  ``bits'' if \code{base = 2}.}

  \item{a}{lower integration limit}

  \item{b}{upper integration limit}
}
\value{
  Entropy (real).
}
\description{
  Computes the Shannon entropy \eqn{H(p)} for a continuous
  probability density function (pdf) \eqn{p(x)}.
}
\details{
  The Shannon entropy of a continuous RV \eqn{X \sim p(x)}
  is defined as \deqn{ H(p) = -\int_{-\infty}^{\infty} p(x)
  \log p(x) d x }

  Contrary to the entropy of a discrete RV, continuous RVs
  can have negative entropy (see Examples).
}
\examples{
my_density <- function(x) {
    dunif(x, -1, 1)
}

\dontrun{
# an error since limits of integration were not specified
continuous_entropy(my_density)
}
continuous_entropy(my_density, -1, 1)  # = log(b - a)

# entropy of the U(a,b) = log(b-a). thus not necessarily positive anymore, e.g.
continuous_entropy(function(x) dunif(x, 0, 0.5), 0, 0.5)  # log2(0.5)
}
\seealso{
  \code{\link{discrete_entropy}}
}
\keyword{math}
\keyword{univar}
