#'@title Estimates spectrum of multivariate time series and returns a 3D array
#'@rdname mvspectrum
#'@aliases mvspectrum normalize_mvspectrum mvspec2mvspectrum SDF2mvspectrum
#'@description
#'Spectrum of a multivariate time series are matrix-valued functions 
#'of the frequency \eqn{\lambda}. \code{mvspectrum} estimates these spectra 
#'and puts them in a 3D array of dimension 
#'\eqn{nfreqs \times nseries \times nseries}.
#'
#'@param series univariate or multivariate time series
#'@param method method for spectrum estimation; see \code{method} argument in
#'\code{\link[sapa]{SDF}} or use \code{'mvspec'} for estimation via 
#'\code{\link[astsa]{mvspec}} from the \pkg{astsa} package.
#'@param normalize indicator; if \code{TRUE} the spectrum will be normalized 
#'such that it adds up to \eqn{0.5} (for univariate spectra) or to 
#'\eqn{0.5 \times I_n} (for multivariate spectra from \eqn{n} time series).
#'@param smoothing indicator; if \code{TRUE} the spectrum gets additionally
#'smoothed using a nonparametric smoother from \code{\link[mgcv]{gam}} with an
#'automatically chosen (cross-validation) smoothing parameter.
#'@param \dots additional arguments passed to \code{\link[sapa]{SDF}} or 
#'\code{\link[astsa]{mvspec}} (e.g. \code{taper})
#'@return 
#'\code{mvspectrum} returns a 3D array of dimension 
#'\eqn{nfreqs \times nseries \times nseries}, where
#'\itemize{
#'  \item nfreqs is the number of frequencies
#'  \item nseries are number of series (number of columns in \code{series}).
#'}
#'@references 
#'See References in \code{\link[stats]{spectrum}}, \code{\link[sapa]{SDF}}, 
#'\code{\link[astsa]{mvspec}}.
#'@keywords ts
#'@export
#'@examples
#'
#'# use SDF first and then SDF2mvspectrum
#'set.seed(1)
#'XX = cbind(rnorm(100), arima.sim(n = 100, list(ar= 0.9)))
#'ss3d = mvspectrum(XX)
#'dim(ss3d)
#'
#'ss3d[2,,] # spectrum at omega_1; in general complex-valued, but Hermitian
#'identical(ss3d[2,,], Conj(t(ss3d[2,,]))) # is Hermitian
#'
mvspectrum <- function(series, method = "multitaper", normalize = FALSE, 
                       smoothing = FALSE, Sigma_X = NULL, ...){
  series = as.matrix(series)
  nseries = ncol(series)
  if (method == "mvspec"){
    out = mvspec2mvspectrum(astsa::mvspec(series, plot = FALSE, ...))
  } else if (method == "ar"){
    out <- spec.ar(c(series), method = "mle", plot = FALSE)$spec[-1]
  } else {
    out = SDF2mvspectrum(sdf.output = SDF(series, method = method, ...))
  }
  
  if (smoothing & nseries == 1){
    mod_temp = gam(log(c(out))~s(c(1:length(out))))
    out = exp(mod_temp$fit)
  }
  
  out = out / 2
  
  if (normalize){
    out = normalize_mvspectrum(out, Sigma_X = Sigma_X)
  }
  
  if (nseries == 1){
    attr(out, "frequency") = seq(0, pi, length = length(out))
  } else {
    attr(out, "frequency") = seq(0, pi, length = nrow(out))
  }
  
  class(out) = "mvspectrum"
  
  invisible(out)
}

#' @rdname mvspectrum
#' @keywords manip
#' @param mvspectrum.output an object of class \code{"mvspectrum"}
#' @param Sigma_X optional; covariance matrix of \code{series}
#' @export
#' @return
#' Returns a normalized spectrum. 
#' \describe{
#'   \item{univariate:}{output adds up to \eqn{0.5}.}
#'   \item{multivariate:}{output adds up to \eqn{0.5 \, I_n}, 
#'         where \eqn{I_n} is the \eqn{n \times n} identitiy matrix.}
#' }
#' 
#' @examples
#' 
#' xx = rnorm(1000)
#' var(xx)
#' mean(mvspectrum(xx, normalize = FALSE, method = "direct")) * 2
#' 
#' mean(mvspectrum(xx, normalize = FALSE, method = "wosa")) * 2
#' mean(mvspectrum(xx, normalize = TRUE, method = "wosa")) * 2
#' 

normalize_mvspectrum <- function(mvspectrum.output, Sigma_X = NULL) {
  nseries <- dim(mvspectrum.output)[2]
  if (is.null(nseries)) {
    nseries <- 1
  }
  
  if (nseries == 1) {
    f3D <- mvspectrum.output/mean(mvspectrum.output)
  } else {
    spectra_sum_inv <- solve(apply(mvspectrum.output, 2:3, mean))
    if (!is.null(Sigma_X)){
      spectra_sum_inv <- spectra_sum_inv %*% Sigma_X
    }
    f3D <- array(t(apply(mvspectrum.output, 1, 
                         function(x) as.matrix(spectra_sum_inv %*% x, ncol = nseries))), 
                 dim(mvspectrum.output))
  }
  invisible(f3D/2)
}
