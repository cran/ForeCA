\name{foreca.EM-aux}
\alias{foreca.EM-aux}
\alias{foreca.EM.E_step}
\alias{foreca.EM.M_step}
\alias{foreca.EM.h}
\title{ForeCA EM auxiliary functions}
\usage{
foreca.EM.E_step(f.U, weightvector)

foreca.EM.M_step(f.U, f.current, base = NULL, minimize = TRUE,
  Sigma.X = NULL)

foreca.EM.h(weightvector.new, f.U, weightvector.current = weightvector.new,
  f.current = NULL, base = NULL)
}
\arguments{
  \item{f.U}{an object of class \code{"mvspectrum"} with
  \code{normalized = TRUE}}

  \item{weightvector}{a weight vector for
  \eqn{y_t=\mathbf{w}' \mathbf{X}_t}. Must have unit norm
  in \eqn{L^2}.}

  \item{base}{logarithm base; if \code{NULL} it sets it
  automatically to \eqn{T}, such that the resulting
  discrete entropy estimate is bounded above by \eqn{1} (it
  is always bounded below by \eqn{0}).}

  \item{minimize}{logical; if \code{TRUE} it returns the
  eigenvector corresponding to the smallest eigenvalue;
  otherwise to the largest eigenvalue.}

  \item{Sigma.X}{optional; covariance matrix of
  \eqn{\mathbf{X}}}

  \item{weightvector.new}{weightvector
  \eqn{\widehat{\mathbf{w}}_{i+1}} of the new iteration
  (i+1)}

  \item{f.current}{spectral density of \eqn{y_t=\mathbf{w}'
  \mathbf{X}_t} for the current estimate
  \eqn{\widehat{\mathbf{w}}_i} (required for
  \code{foreca.EM.M_step}; optional for
  \code{foreca.EM.h}).}

  \item{weightvector.current}{weightvector
  \eqn{\widehat{\mathbf{w}}_{i}} of the current iteration
  (i)}
}
\value{
\code{foreca.EM.E_step} returns a vector containing the
normalized spectral density (normalized such that its
\code{mean} equals \eqn{0.5} - since \code{f.U} only
contains positive frequencies and is symmetric).

\code{foreca.EM.M_step} returns a list with three elements:
\describe{ \item{\code{matrix}:}{is the weighted covariance
matrix (positive semi-definite), where the weights are the
negative log of the spectral density,}
\item{\code{vector}:}{minimizing (or maximizing if
\code{minimize = FALSE}) eigenvector of \code{matrix},}
\item{\code{value}:}{corresponding eigenvalue.} }

\code{foreca.EM.h} returns (see References for details):
\itemize{ \item the negative entropy (non-negative real) if
\code{weightvector.new = weightvector.current} \item an
upper bound of that entropy for \code{weightvector.new} if
\code{weightvector.new != weightvector.current} }
}
\description{
The EM-like algorithm relies on several auxiliary
functions:

\code{foreca.EM.E_step} computes the spectral density of
\eqn{y_t=\mathbf{w}' \mathbf{X}_t} given the weightvector
\eqn{\mathbf{w}}.

\code{foreca.EM.M_step} computes the minimizing eigenvector
(\eqn{\rightarrow \widehat{\mathbf{w}}_{i+1}}) of the
weighted covariance matrix, where the weights equal the
negative logarithm of the spectral density at the current
\eqn{\widehat{\mathbf{w}}_i}.

\code{foreca.EM.h} evaluates the entropy of the spectral
density as a function of \eqn{\mathbf{w}}. This is the
objective funcion that should be minimized.
}
\examples{
XX <- diff(log(EuStockMarkets)) * 100
ff <- mvspectrum(XX, 'wosa', normalize = TRUE)

ww0 <- matrix(rnorm(ncol(XX)))
ww0 <- ww0 / norm(ww0, 'F')

f.ww0 <- foreca.EM.E_step(ff, ww0)
onestep <- foreca.EM.M_step(ff, f.ww0)
image(onestep$matrix)
\dontrun{
require(LICORS)
# if you have 'LICORS' package installed, better use this:
image2(onestep$matrix)
}
ww1 <- onestep$vector
f.ww1 <- foreca.EM.E_step(ff, ww1)

layout(matrix(1:2, ncol = 2))
#par(mar = c(4, 4, 1, 1), mfcol= c(1,2))
matplot(seq(0, pi, length = length(f.ww0)), cbind(f.ww0, f.ww1),
        type = "l", lwd =2, xlab = "omega_j", ylab = "f(omega_j)")
plot(f.ww0, f.ww1, pch = ".", cex = 3, xlab = "iteration 0",
     ylab = "iteration 1", main = "Spectral density")
abline(0, 1, col = 'blue', lty = 2)
foreca.EM.h(ww0, ff)       # iteration 0
foreca.EM.h(ww1, ff, ww0)  # min eigenvalue inequality
foreca.EM.h(ww1, ff)       # KL divergence inequality
onestep$value

Omega(spectrum.estimate = f.ww0) / 100 + foreca.EM.h(ww0, ff)
Omega(spectrum.estimate = f.ww1) / 100 + foreca.EM.h(ww1, ff)
}
\keyword{manip}

