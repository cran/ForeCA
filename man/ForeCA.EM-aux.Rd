\name{foreca.EM-aux}
\alias{foreca.EM-aux}
\alias{foreca.EM.E_step}
\alias{foreca.EM.h}
\alias{foreca.EM.init_weightvector}
\alias{foreca.EM.M_step}
\title{ForeCA EM auxiliary functions}
\usage{
foreca.EM.E_step(f_U, weights = NULL)

foreca.EM.M_step(f_U, f_current = NULL, base = NULL, minimize = TRUE, Sigma_X = NULL)

foreca.EM.h(weights_new, f_U, weights_current = weights_new, f_current = NULL, base = NULL)

foreca.EM.init_weightvector(UU, f_U, method = "SFA", seed = NULL, lag = 1)
}
\arguments{
  \item{f_U}{an object of class \code{"mvspectrum"} with
  \code{normalized = TRUE}}

  \item{weights}{a weight vector for \eqn{y_t=\mathbf{w}'
  \mathbf{X}_t}. Must have unit norm in \eqn{L_2}.}

  \item{base}{logarithm base; if \code{NULL} it sets it
  automatically to \eqn{T}, such that the resulting entropy
  estimate is bounded above by \eqn{1} (it is always
  bounded below by \eqn{0}).}

  \item{minimize}{indicator; if \code{TRUE} it returns the
  eigenvector correpsonding to the smallest eigenvalue;
  otherwise to the largest eigenvalue.}

  \item{Sigma_X}{optional; covariance matrix of
  \eqn{\mathbf{X}}}

  \item{weights_new}{weights
  \eqn{\widehat{\mathbf{w}}_{i+1}} of the new iteration
  (i+1)}

  \item{f_current}{spectral density of the current estimate
  \eqn{\widehat{\mathbf{w}}_i} (required for
  \code{foreca.EM.M_step}; optional for
  \code{foreca.EM.h}).}

  \item{weights_current}{weights
  \eqn{\widehat{\mathbf{w}}_{i}} of the current iteration
  (i)}

  \item{UU}{uncorrelated multivariate time series}

  \item{method}{string indicating the method to estimate
  the initial weightvector; default \code{"SFA"}. See
  Details.}

  \item{seed}{seed to use for random initialization; if
  \code{NULL}, it sets a random seed}

  \item{lag}{integer; lag for the difference operator;
  default \code{lag=1}.}
}
\value{
  \code{foreca.EM.E_step} returns a vector containing the
  normalized spectral density (normalized such that its
  \code{mean} equals \eqn{0.5} - since \code{f_U} only
  contains positive frequencies and is symmetric).

  \code{foreca.EM.M_step} returns a list with three
  elements: \describe{ \item{\code{matrix}:}{is the
  weighted covariance matrix (positive semi-definite),
  where the weights are the negative log of the spectral
  density} \item{\code{vector}:}{minimizing (or maximizing
  if \code{minimize = FALSE}) eigenvector of that matrix}
  \item{\code{value}:}{corresponding eigenvalue} }

  \code{foreca.EM.h} returns (see References for details):
  \itemize{ \item the negative entropy (non-negative real)
  if \code{weights_new = weights_current} \item an upper
  bound of that entropy for \code{weights_new} if
  \code{weights_new != weights_current} }

  \code{foreca.EM.init_weightvector} returns a good
  starting vector \eqn{\mathbf{w}_0} for
  \code{\link{foreca.EM.opt_comp}}.
}
\description{
  The EM-like algorithm needs three auxiliary functions:

  \code{foreca.EM.E_step} computes the spectral density
  given a weight-vector \eqn{\mathbf{w}_i}.

  \code{foreca.EM.M_step} computes the minimizing
  eigenvector (\eqn{\rightarrow
  \widehat{\mathbf{w}}_{i+1}}) of the weighted covariance
  matrix, where the weights equal the negative logarithm of
  the spectral density at the current
  \eqn{\widehat{\mathbf{w}}_i}.

  \code{foreca.EM.h} evaluates the entropy of the spectral
  density as a function of \eqn{\mathbf{w}}. This is the
  objective funcion that should be minimized.

  \code{foreca.EM.init_weightvector} returns an
  initialization vector \eqn{\mathbf{w}_0 \in R^k}, with
  \eqn{L_2} norm \eqn{1}, for the EM-like algorithm in
  \code{\link{foreca.EM.opt_comp}}. They are all based on
  quickly computable heuristics than approximate a
  forecastable signal. See \code{method} argument and
  Details below.
}
\details{
  The \code{method} argument determines what heuristics
  should be used to get a good starting vector
  \eqn{\mathbf{w}_0}. This vector is length \eqn{k} and has
  unit norm in \eqn{L_2}:

  \describe{ \item{max}{uses the vector with all \eqn{0}s,
  but a \eqn{1} at the position of the maximum
  forecastability of the series in \code{UU}.}
  \item{SFA_slow}{uses the first eigenvector of SFA
  (slowest signal).} \item{SFA_fast}{uses the last
  eigenvector of SFA (fastes signal).} \item{SFA}{checks
  both first and last, and chooses the one with higher
  forecastability.} \item{cauchy}{random starts using
  \code{rcauchy(k)}} \item{unif}{random starts using
  \code{runif(k, -1, 1)}} \item{norm}{random starts using
  \code{runif(k, 0, 1)}} }

  Slow Feature Analysis (SFA; see References) finds
  \emph{slow} signals. For details see the References; for
  us it's important that SFA is equivalent to finding the
  signal with largest lag \eqn{1} autocorrelation.  This is
  not necessarily the most forecastable, but a good start.

  The disadvantage of SFA for forecasting is that e.g.
  white noise (WN) is ranked higher than an AR(1) with
  negative autocorrelation.  While it is true that WN is
  slower, it is not more forecastable.  Thus we are also
  interested in the fastest signal, i.e. the last
  eigenvector.  determine a good starting vector. The
  slowest feature corresponds to maximizing the lag \eqn{1}
  correlation
}
\examples{
XX <- diff(log(EuStockMarkets)) * 100
ff <- mvspectrum(XX, "wosa", normalize = TRUE)

ww0 <- matrix(rnorm(ncol(XX)))
ww0 <- ww0/norm(ww0, "F")

f_ww0 <- foreca.EM.E_step(ff, ww0)
onestep <- foreca.EM.M_step(ff, f_ww0)
image(onestep$matrix)
\dontrun{
# if you have 'LICORS' package installed, better use this:
image2(onestep$matrix)
}
ww1 <- onestep$vector
f_ww1 <- foreca.EM.E_step(ff, ww1)

par(mar = c(4, 4, 1, 1), mfcol = c(1, 2))
matplot(seq(0, pi, length = length(f_ww0)), cbind(f_ww0, f_ww1), type = "l", lwd = 2, xlab = "frequencies")
plot(f_ww0, f_ww1, pch = ".", cex = 3, xlab = "iteration 0", ylab = "iteration 1")
foreca.EM.h(ww0, ff)  # iteration 0
foreca.EM.h(ww1, ff, ww0)  # min eigenvalue inequality
foreca.EM.h(ww1, ff)  # KL divergence inequality
onestep$value

Omega(spectrum_estimate = f_ww0)/100 + foreca.EM.h(ww0, ff)
Omega(spectrum_estimate = f_ww1)/100 + foreca.EM.h(ww1, ff)
}
\references{
  Laurenz Wiskott and Terrence J. Sejnowski (2002).
  \dQuote{Slow Feature Analysis: Unsupervised Learning of
  Invariances}, Neural Computation 14:4, 715-770.
}
\keyword{manip}
