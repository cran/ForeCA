#'@rdname mvspectrum
#'
#'@description
#'\code{mvspec2mvspectrum} converts output from the \code{\link[astsa]{mvspec}} 
#'function to the \code{\link{mvspectrum}} output.
#'@details
#'The \code{\link[astsa]{mvspec}} function has frequencies in the first dimension,
#'\code{\link{mvspectrum}} in the last.  \code{mvspec2mvspectrum} simply 
#'reshapes this \eqn{3D} array to an array of size 
#'\eqn{nfreqs \times nseries \times nseries}.
#'
#'@param mvspec.output output from \code{\link[astsa]{mvspec}}
#'@return 
#'An array of dimension \eqn{nfreqs \times nseries \times nseries}.
#'@keywords ts manip

mvspec2mvspectrum = function(mvspec.output){
  out = R.utils::wrap(mvspec.output$fxx, map = list(3,1,2))
  attr(out, "frequencies") = mvspec.output$freq
  attr(out, "spectrum") = mvspec.output$spec
  invisible(out)
}


