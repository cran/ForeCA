\name{foreca-utils}
\alias{biplot.foreca}
\alias{foreca-utils}
\alias{plot.foreca}
\alias{plot.foreca.EM.opt_weightvector}
\alias{summary.foreca}
\title{Plot, summary, and print methods for class 'foreca'}
\usage{
\method{plot}{foreca}(x, lag = 10, alpha = 0.05, ...)

\method{plot}{foreca.EM.opt_weightvector}(x, main = "", cex.lab = 1.1, ...)

\method{summary}{foreca}(object, lag = 10, alpha = 0.05, ...)

\method{biplot}{foreca}(x, ...)
}
\arguments{
  \item{...}{additional arguments passed to
  \code{\link[stats]{biplot.princomp}},
  \code{\link[stats]{biplot.default}},
  \code{\link[graphics]{plot}}, or
  \code{\link[base]{summary}}.}

  \item{x}{an object of class \code{"foreca"} or
  \code{"foreca.EM.opt_weightvector"}}

  \item{main}{an overall title for the plot: see
  \code{\link[graphics]{title}}.}

  \item{cex.lab}{size of the axes labels}

  \item{object}{an object of class \code{"foreca"}}

  \item{alpha}{significance level for testing white noise
  in \code{\link[stats]{Box.test}}; default: \eqn{5\%}.}

  \item{lag}{how many lags to test in
  \code{\link[stats]{Box.test}}; default \eqn{10} lags.}
}
\description{
Plot, summary, and print methods for class 'foreca'

\code{plot.foreca} shows a visual summary of the ForeCA
results with biplots, screeplots, and white noise tests.

\code{plot.foreca.EM.opt_weightvector} shows how the
EM-like algorithm arrived at the i-th optimal a
weight-vector \eqn{\mathbf{w}_i^*}.  It shows trace plots
of the objective function (\code{\link{foreca.EM.h}}) and
of the weight vector, and the transformed signal
\eqn{y_t^*} along with its spectral density estimate
\eqn{\widehat{f}_y(\omega_j)}.

\code{summary.foreca} computes summary statistics of the
ForeCA results.

\code{biplot.foreca} shows a biplot of the ForeCA
weightvectors (a wrapper around
\code{\link[stats]{biplot.princomp}}).
}
\examples{
# see examples in 'foreca.EM'
}
\keyword{hplot}
\keyword{manip}

