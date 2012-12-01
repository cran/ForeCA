#'@rdname mvspectrum
#'
#'@description
#'\code{SDF2mvspectrum} converts \code{\link[sapa]{SDF}} output to a 3D array with
#'number of frequencies in the first dimension and the spectral density matrix
#'in the other two (see Details below).
#'
#'@details
#'The \code{\link[sapa]{SDF}} function only returns the upper
#'diagonal (including diagonal), since spectrum matrices are Hermitian. For
#'efficient vectorized computations, however, the full matrices are required.
#'Thus \code{SDF2mvspectrum} converts SDF output to a 3D array with
#'number of frequencies in the first dimension and the spectral density matrix
#'in the other two (see Details below).
#'
#'\code{SDF2mvspectrum} is typically not called by the user, but by
#'\code{\link{mvspectrum}}.
#'
#'@param sdf.output an object of class \code{"SDF"} from \code{\link[sapa]{SDF}}
#'@keywords ts manip

SDF2mvspectrum <-
function(sdf.output){
  nseries = attr(sdf.output, "n.series")
  nobs = attr(sdf.output, "n.sample")
  nfreqs = length(attr(sdf.output, "frequency"))
  f_lambda = array(0, dim = c(nfreqs, nseries, nseries))
  
  if (nseries > 1){

    col.sels = 1:nseries
    ntotal = nseries
    
    for (ii in 1:ntotal){
      f_lambda[, ii, ii:ntotal] = sdf.output[, col.sels]
      col.sels = col.sels + (ntotal - ii)
      col.sels = col.sels[-1]
    }
    
    fill_symmetric = function(mat){
      mat = mat + Conj(t(mat))
      ind <- lower.tri(mat)
      mat[ind] <- t(Conj(mat))[ind] 
      return(mat)
    }
      
    f_lambda = apply(f_lambda, 1, fill_symmetric)
    f_lambda = array(t(f_lambda), dim = c(nfreqs, nseries, nseries))
  } else {
    f_lambda[, 1, 1] = c(sdf.output)
  }
  f_lambda = f_lambda[-1, , ]
  attr(f_lambda, "frequency") = attr(sdf.output, "frequency")[-1]

  invisible(f_lambda)
}
