\name{ForeCA-utils}
\alias{biplot.ForeCA}
\alias{ForeCA-utils}
\alias{plot.ForeCA}
\alias{plot.ForeCA.EM.one_weightvector}
\alias{summary.ForeCA}
\title{ForeCA plots, summary, print, etc methods}
\usage{
\method{plot}{ForeCA}(x, lag = 10, alpha = 0.05, ...)

\method{plot}{ForeCA.EM.one_weightvector}(x, main = "", ...)

\method{summary}{ForeCA}(object, lag = 10, alpha = 0.05, ...)

\method{biplot}{ForeCA}(x, ...)
}
\arguments{
  \item{...}{additional arguments passed to
  \code{\link[stats]{biplot.princomp}},
  \code{\link[stats]{biplot.default}},
  \code{\link[graphics]{plot}}, or
  \code{\link[base]{summary}}.}

  \item{x}{an object of class \code{"ForeCA"} or
  \code{"ForeCA.EM.one_weightvector"}}

  \item{main}{an overall title for the plot: see
  \code{\link[graphics]{title}}.}

  \item{object}{an object of class \code{"ForeCA"}}

  \item{alpha}{significance level for testing white noise
  in \code{\link[stats]{Box.test}}}

  \item{lag}{how many lags to test in
  \code{\link[stats]{Box.test}}}
}
\description{
  ForeCA plots, summary, print, etc methods

  \code{plot.ForeCA} plots the output of a \code{"ForeCA"}
  object. It is a visual summary of the analysis with
  biplots, screeplots, and white noise tests.

  \code{plot.ForeCA.EM.one_weightvector} gives a visual
  summary of the EM-like algorithm to obtain a
  weight-vector \eqn{\mathbf{w}_i^*}.  It gives trace plots
  of the objective function (\code{\link{ForeCA.EM.h}}) and
  of the weight vector, and it shows the transformed signal
  \eqn{y_t^*} along with its spectral density estimate
  \eqn{\widehat{f}_y(\omega_j)}.

  \code{summary.ForeCA} is a \code{\link[base]{summary}}
  method for a \code{\link{ForeCA}} output.

  \code{biplot.ForeCA} plots a biplot from the output of
  \code{\link{ForeCA}} (a wrapper around
  \code{\link[stats]{biplot.princomp}}).
}
\examples{
# see examples in 'ForeCA.EM'
}
\keyword{hplot}
\keyword{manip}
