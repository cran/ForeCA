\name{foreca-utils}
\alias{biplot.foreca}
\alias{foreca-utils}
\alias{plot.foreca}
\alias{plot.foreca.EM.opt_comp}
\alias{summary.foreca}
\title{Plots, summary, print, etc methods for class 'foreca'}
\usage{
\method{plot}{foreca}(x, lag = 10, alpha = 0.05, ...)

\method{plot}{foreca.EM.opt_comp}(x, main = "", ...)

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
  \code{"foreca.EM.opt_comp"}}

  \item{main}{an overall title for the plot: see
  \code{\link[graphics]{title}}.}

  \item{object}{an object of class \code{"foreca"}}

  \item{alpha}{significance level for testing white noise
  in \code{\link[stats]{Box.test}}}

  \item{lag}{how many lags to test in
  \code{\link[stats]{Box.test}}}
}
\description{
  Plots, summary, print, etc methods for class 'foreca'

  \code{plot.foreca} plots the output of a \code{"foreca"}
  object. It is a visual summary of the analysis with
  biplots, screeplots, and white noise tests.

  \code{plot.foreca.EM.opt_comp} gives a visual summary of
  the EM-like algorithm to obtain a weight-vector
  \eqn{\mathbf{w}_i^*}.  It gives trace plots of the
  objective function (\code{\link{foreca.EM.h}}) and of the
  weight vector, and it shows the transformed signal
  \eqn{y_t^*} along with its spectral density estimate
  \eqn{\widehat{f}_y(\omega_j)}.

  \code{summary.foreca} is a \code{\link[base]{summary}}
  method for an object of class \code{"foreca"}.

  \code{biplot.foreca} plots a biplot for an object of
  class \code{"foreca"} (a wrapper around
  \code{\link[stats]{biplot.princomp}}).
}
\examples{
# see examples in 'foreca.EM'
}
\keyword{hplot}
\keyword{manip}
