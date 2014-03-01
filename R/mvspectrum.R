#'@title Estimates spectrum of multivariate time series
#'@rdname mvspectrum
#'@aliases mvspectrum normalize_mvspectrum mvspec2mvspectrum SDF2mvspectrum
#'@description
#'Spectra of a multivariate time series are matrix-valued functions 
#'of the frequency \eqn{\lambda}. \code{mvspectrum} estimates these spectra 
#'and puts them in a 3D array of dimension \eqn{num.freqs \times K \times K}.
#'
#'@param series univariate or multivariate time series
#'@param method method for spectrum estimation; see \code{method} argument in
#'\code{\link[sapa]{SDF}} or use \code{'mvspec'} for estimation via 
#'\code{\link[astsa]{mvspec}} from the \pkg{astsa} package.
#'@param normalize logical; if \code{TRUE} the spectrum will be normalized 
#'such that it adds up to \eqn{0.5} (for univariate spectra) or to 
#'\eqn{0.5 \times I_n} (for multivariate spectra from \eqn{n} time series).
#'@param smoothing logical; if \code{TRUE} the spectrum gets additionally
#'smoothed using a nonparametric smoother from \code{\link[mgcv]{gam}} with an
#'automatically chosen (cross-validation) smoothing parameter.
#'@param \dots additional arguments passed to \code{\link[sapa]{SDF}} or 
#'\code{\link[astsa]{mvspec}} (e.g., \code{taper})
#'@return 
#'\code{mvspectrum} returns a 3D array of dimension 
#'\eqn{num.freqs \times K \times K}, where
#'\itemize{
#'  \item num.freqs is the number of frequencies
#'  \item K is the number of series (columns in \code{series}).
#'}
#'@references 
#'See References in \code{\link[stats]{spectrum}}, \code{\link[sapa]{SDF}}, 
#'\code{\link[astsa]{mvspec}}.
#'@keywords ts
#'@export
#'@importFrom sapa SDF
#'@importFrom mgcv gam
#'@importFrom ifultools checkScalarType
#'@import nlme
#'@import splus2R
#'@examples
#'
#'# use SDF first and then SDF2mvspectrum
#'set.seed(1)
#'XX <- cbind(rnorm(100), arima.sim(n = 100, list(ar= 0.9)))
#'ss3d <- mvspectrum(XX)
#'dim(ss3d)
#'
#'ss3d[2,,] # at omega_1; in general complex-valued, but Hermitian
#'identical(ss3d[2,,], Conj(t(ss3d[2,,]))) # is Hermitian
#'

mvspectrum <- function(series, 
                       method = 
                         c("multitaper", "direct", "lag window", "wosa", 
                           "mvspec", "ar"),
                       normalize = FALSE, smoothing = FALSE, Sigma.X = NULL, ...) {
  
  method <- match.arg(method)
  series <- as.matrix(series)
  num.series <- ncol(series)
  num.obs <- nrow(series)
  
  if (method == "mvspec") {
    out <- mvspec2mvspectrum(astsa::mvspec(series, plot = FALSE, 
                                           detrend = FALSE, fast = FALSE, ...))
  } else if (method == "ar") {
    stopifnot(num.series == 1)
    out <- spec.ar(c(series), method = "burg", plot = FALSE,
                   n.freq = ceiling(length(series) / 2 + 1))$spec[-1]
  } else {
    out <- SDF2mvspectrum(sdf.output = SDF(series, method = method, 
                                           recenter = TRUE, ...))
    out <- as.array(out)
    if (length(dim(out)) == 1) {
      out <- array(out, dim = c(length(out), 1, 1))
    }
    # keep only every second frequency starting at first non-zero frequency
    # (SDF returns spectral densities on [0, pi] of length T; and not of T/2)
    every.2nd <- seq(2, dim(out)[1], by = 2)
    out <- out[every.2nd,,]
  }
  
  if (smoothing && num.series == 1) {
    mod.tmp <- gam(c(out) ~ s(c(seq_along(out))),
                   family = gaussian(), method = "REML")
    out <- mod.tmp$fit
  }
  
  out <- out/2
  
  if (normalize) {
    out <- normalize_mvspectrum(out, Sigma.X = Sigma.X)
  }
  # add frequencies from (0, pi]; remove '0' frequency
  if (num.series == 1) {
    num.freqs.in.0.pi <- length(out)
  } else {
    num.freqs.in.0.pi <- nrow(out)
  }
  attr(out, "frequency") <- 
    seq(1, num.freqs.in.0.pi, by = 1) / (2 * num.freqs.in.0.pi + 1) * (2 * pi)
  class(out) <- "mvspectrum"
  
  invisible(out)
} 
#' @rdname mvspectrum
#' @keywords manip
#' @param mvspectrum.output an object of class \code{"mvspectrum"}
#' @param Sigma.X optional; covariance matrix of \code{series}
#' @export
#' @return
#' \code{normalize_mvspectrum} returns a normalized spectrum:
#' \describe{
#'   \item{univariate:}{output adds up to \eqn{0.5}.}
#'   \item{multivariate:}{output adds up to \eqn{0.5 \, I_n}, 
#'         where \eqn{I_n} is the \eqn{n \times n} identity matrix.}
#' }
#' 
#' @examples
#' 
#' xx <- rnorm(1000)
#' var(xx)
#' mean(mvspectrum(xx, normalize = FALSE, method = "direct")) * 2
#' 
#' mean(mvspectrum(xx, normalize = FALSE, method = "wosa")) * 2
#' mean(mvspectrum(xx, normalize = TRUE, method = "wosa")) * 2
#' 

normalize_mvspectrum <- function(mvspectrum.output, Sigma.X = NULL) {
  num.series <- dim(mvspectrum.output)[2]
  if (is.null(num.series)) {
    num.series <- 1
  }
  
  if (num.series == 1) {
    f3D <- mvspectrum.output/mean(mvspectrum.output)
  } else {
    spectra.sum.inv <- solve(apply(mvspectrum.output, 2:3, mean))
    if (!is.null(Sigma.X)) {
      spectra.sum.inv <- spectra.sum.inv %*% Sigma.X
    }
    f3D <- array(t(apply(mvspectrum.output, 1, 
                         function(x) as.matrix(spectra.sum.inv %*% 
                                                 x, ncol = num.series))), 
                 dim(mvspectrum.output))
  }
  invisible(f3D / 2)
} 