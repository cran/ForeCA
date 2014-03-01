#'@rdname mvspectrum
#'
#'@description
#'\code{mvspec2mvspectrum} converts output from the \code{\link[astsa]{mvspec}} 
#'function to the \code{\link{mvspectrum}} output.
#'@details
#'The \code{\link[astsa]{mvspec}} function returns the multivariate spectrum
#'in a \eqn{3D} array with frequencies in the first dimension, whereas 
#'\code{\link{mvspectrum}} in the last.  \code{mvspec2mvspectrum} simply 
#'reshapes the former to the latter array.
#'
#'@param mvspec.output output from \code{\link[astsa]{mvspec}}
#'@return 
#'\code{mvspec2mvspectrum} returns the same object as \code{\link{mvspectrum}}.
#'@keywords ts manip
#'@export

mvspec2mvspectrum <- function(mvspec.output) {
  out <- R.utils::wrap(mvspec.output$fxx, map = list(3, 1, 2))
  if (ncol(out) == 1) {
    # remove '+0i' imaginary part for univariate spectra
    out <- Re(out)
  }
  attr(out, "frequencies") <- mvspec.output$freq
  attr(out, "spectrum") <- mvspec.output$spec
  invisible(out)
} 

