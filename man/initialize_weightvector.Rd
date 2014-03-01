\name{initialize_weightvector}
\alias{initialize_weightvector}
\title{Initialize the weightvector for iterative ForeCA algorithms}
\usage{
initialize_weightvector(U, f.U, method = c("SFA", "SFA.slow", "SFA.fast",
  "cauchy", "unif", "norm", "max"), seed = NULL, lag = 1)
}
\arguments{
  \item{f.U}{an object of class \code{"mvspectrum"} with
  \code{normalized = TRUE}}

  \item{U}{uncorrelated multivariate time series of
  dimension \eqn{T \times K}.}

  \item{method}{string indicating the method to estimate
  the initial weightvector; default \code{"SFA"}; see
  Details.}

  \item{seed}{seed to use for random initialization; if
  \code{NULL}, it sets a random seed (and it will be
  returned for reproducibility).}

  \item{lag}{integer; lag for the difference operator;
  default \code{lag=1}.}
}
\value{
A vector of length \eqn{k} with unit norm (in \eqn{L^2}).
}
\description{
\code{initialize_weightvector} returns a unit norm (in
\eqn{L^2}) vector \eqn{\mathbf{w}_0 \in R^k} that can be
used as the initial starting point for any iterative ForeCA
algorithm, e.g., the EM-like algorithm in
\code{\link{foreca.EM.opt_weightvector}}. Several quickly
computable heuristics to get a good starting value are
available. See \code{method} argument and Details below.
}
\details{
The \code{method} argument specifies the heuristics that is
used to get a good starting vector \eqn{\mathbf{w}_0}. This
vector has length \eqn{k} and unit norm in \eqn{L^2}:

\describe{ \item{max}{uses the vector with all \eqn{0}s,
but a \eqn{1} at the position of the maximum forecastable
series in \code{U}.} \item{SFA.slow}{uses the first
eigenvector of SFA (slowest signal).} \item{SFA.fast}{uses
the last eigenvector of SFA (fastest signal).}
\item{SFA}{checks both slow and fast, and chooses the one
with higher forecastability.} \item{cauchy}{random starts
using \code{rcauchy(k)}} \item{unif}{random starts using
\code{runif(k, -1, 1)}} \item{norm}{random starts using
\code{rnorm(k, 0, 1)}} }

Slow Feature Analysis (SFA) finds \emph{slow} signals (see
References below), and can be quickly (and analytically)
computed solving a generalized eigen-value problem.  For
ForeCA it is important to know that SFA is equivalent to
finding the signal with largest lag \eqn{1}
autocorrelation.  This is not necessarily the most
forecastable, but a good start.

The disadvantage of SFA for forecasting is that e.g., white
noise (WN) is ranked higher than an AR(1) with negative
autocorrelation coefficient \eqn{\rho_1 < 0}.  While it is
true that WN is slower, it is not more forecastable.  Thus
we are also interested in the fastest signal, i.e., the
last eigenvector. The so obtained fastest signal
corresponds to minimizing the lag 1 auto-correlation
(possibly \eqn{\rho_1 < 0}).
}
\references{
Laurenz Wiskott and Terrence J. Sejnowski (2002).
\dQuote{Slow Feature Analysis: Unsupervised Learning of
Invariances}, Neural Computation 14:4, 715-770.
}
\keyword{manip}

