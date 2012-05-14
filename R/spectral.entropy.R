spectral.entropy <-
function# Compute the spectral entropy of a time series
  (series, 
   ### univariate time series; multivariate also supported but not really meaningful
   spectrum_method = "wosa", 
   ### specify the method for the spectrum estimation; 
   ### see \code{method} argument in \code{\link[sapa]{SDF}}
   base = exp(1), 
   ### base of the logarithm
   entropy_method = "MLE", 
   ### method to compute the entropy from a sequence of discrete probabilities;
   ### here those 'probabilities' is the spectral density evaluated at the
   ### Fourier frequencies.
   threshold = 0, 
   ### A threshold to set spectrum to 0 if it does not exceed that threshold
   ...
   ### additional arguments
   ){
  
  series = as.matrix(series)
  nseries = ncol(series)

  if (nseries > 1) {
    sdf = spectrum.3D.array(series, method = spectrum_method, ...)
  } else {
    sdf = SDF(series, method = spectrum_method, ...)
  }
  
  if (nseries > 1){
    nfreqs = nrow(sdf)
    spec.ent = sdf[1,,] * log(sdf[1,,], base = base)
    # improve this with an apply command
    for (ff in 2:nfreqs){
      spec.ent = spec.ent + sdf[ff,,] %*% log(sdf[ff,,], base = base)
    }
    spec.ent = spec.ent + Conj(spec.ent) - 
               sdf[1,,] %*% log(sdf[1,,], base = base)
  } else {
    spec.ent = discrete.entropy(c(rev(c(sdf)), c(sdf)), base = base, 
                                method = entropy_method,
 					                      threshold = threshold)
  }
  return(spec.ent)
  ### a real value for the spectral entropy \eqn{H_s(x_t)}.
}
