\name{discrete_entropy}
\alias{discrete_entropy}
\title{Shannon entropy for discrete pmf}
\usage{
discrete_entropy(probs, base = 2, method = c("MLE", "smoothed"),
  threshold = 0, smoothing.weight = 1)
}
\arguments{
  \item{probs}{probabilities (empirical frequencies). Must
  add up to \eqn{1}.}

  \item{base}{see \code{base} argument of
  \code{\link{continuous_entropy}}}

  \item{method}{method to estimate entropy; see Details
  below.}

  \item{threshold}{frequencies below \code{threshold} are
  set to \eqn{0}; default \code{threshold = 0}.}

  \item{smoothing.weight}{how much weight does the uniform
  smoothing term get? default \code{1}.}
}
\value{
Entropy (a non-negative real value).
}
\description{
Computes the Shannon entropy \eqn{\mathcal{H}(p)} of a
discrete RV \eqn{X} taking values in \eqn{\lbrace x_1,
\ldots, x_n \rbrace} with probability mass function (pmf)
\eqn{Prob(X = x_i) \in \lbrace p_1, \ldots, p_n \rbrace}.
}
\details{
\code{\link{discrete_entropy}} uses a plug-in estimator
(\code{method = ''MLE''}): \deqn{ \widehat{H}(p) = -
\sum_{i=1}^{n} \widehat{p}_i \log \widehat{p}_i. }

The smoothed estimator (\code{method = 'smoothed'}) adds a
uniform distribution \eqn{1/n} to every probability (and
then normalizes again): \deqn{ \widehat{p}_i \leftarrow
\widehat{p_i} + \alpha \cdot \frac{1}{n}, \quad i=1,
\ldots, n, } where \eqn{\alpha} is the smoothing weight
parameter.
}
\examples{
pmfs <- rexp(5)
pmfs <- sort(pmfs / sum(pmfs))

unif.distr <- rep(1/length(pmfs), length(pmfs))

matplot(cbind(pmfs, unif.distr), pch = 19, ylab = "P(X = k)", xlab = "k")
matlines(cbind(pmfs, unif.distr))
legend("topleft", c("non-uniform","uniform"), pch = 19, lty = 1:2, col = 1:2)

discrete_entropy(pmfs)
# uniform has largest entropy among all bounded distributions (here = log(5))
discrete_entropy(unif.distr)
# no uncertainty if one element occurs with probability 1
discrete_entropy(c(1,0,0,0,0))
}
\seealso{
\code{\link{continuous_entropy}}
}
\keyword{math}
\keyword{univar}

