\name{discrete.entropy}
\alias{discrete.entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Entropy estimator for a discrete RV given empirical frequencies }
\description{
Estimates the entropy given empirical frequencies \eqn{\widehat{p}_i, i = 1, \ldots, n}.
}
\usage{
discrete.entropy(prob, base = exp(1), method = "MLE", threshold = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prob}{ probabilities (empirical frequencies) }
  \item{base}{ logarithm base; default \code{base = e}. Entropy is measures in ``nats'' for \code{base = e}; in ``bits'' if \code{base = 2}.}
  \item{method}{ method to estimate entropy; see Details below.}
  \item{threshold}{ frequencies below \code{threshold} are set to \eqn{0}; default \code{threshold = 0}. }
}

\details{
Since we only have the frequencies (which equal the spectral density \eqn{\widehat{f}(\omega_j)}), but not the observed counts, currently \code{discrete.entropy} only uses a plug-in estimator (\code{method = ''MLE''}):
\deqn{
\widehat{H}(p) = - \sum_{i=1}^{n} \widehat{p}_i \log \widehat{p}_i.
}

A smoothed estimator is also available (\code{method = ''smoothed''}), which adds \eqn{1/n} to every empirical frequency (and then normalizes again):
\deqn{
\widehat{p}_i \leftarrow \widehat{p_i} + \frac{1}{n}, \quad i=1, \ldots, n.
}

}

\value{
Entropy for observed frequencies (a non-negative real value).
}

\author{
Georg M. Goerg
}

\examples{
probs = rexp(5)
probs = sort(probs / sum(probs))

unif.distr = rep(1/length(probs), length(probs))

matplot(cbind(probs, unif.distr), pch = 19, ylab = "P(X = k)", xlab = "k")
matlines(cbind(probs, unif.distr))
legend("topleft", c("non-uniform","uniform"), pch = 19, lty = 1:2, col = 1:2)

discrete.entropy(probs)
discrete.entropy(unif.distr) # uniform distribution has largest entropy among all distribution
discrete.entropy(c(1,0,0,0,0)) # no uncertainty if one element happens with prob 1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
\keyword{univar}
