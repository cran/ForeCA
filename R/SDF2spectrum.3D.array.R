SDF2spectrum.3D.array <-
function(sdf.output, frequency_0 = TRUE){
  nseries = attr(sdf.output, "n.series")
  nobs = attr(sdf.output, "n.sample")
  nfreqs = length(attr(sdf.output, "frequency"))
  
  f_lambda = array(0, dim = c(nfreqs, nseries, nseries))
  col.sels = 1:nseries
  ntotal = nseries
  
  for (ii in 1:ntotal){
    f_lambda[, ii, ii:ntotal] = sdf.output[, col.sels]
    col.sels = col.sels + (ntotal - ii)
    col.sels = col.sels[-1]
  }
  
  fill_symmetric = function(mat){
    mat = mat + Conj(t(mat))
    #mat = mat - diag(diag(mat)/2)
    ind <- lower.tri(mat)
    mat[ind] <- t(Conj(mat))[ind] 
    #mat[lower.tri(mat)] = Conj(t(mat[upper.tri(mat)]))
    return(mat)
  }
    
  f_lambda = apply(f_lambda, 1, fill_symmetric)
  f_lambda = array(t(f_lambda), dim = c(nfreqs, nseries, nseries))
   
  if (!frequency_0) {
    f_lambda = f_lambda[-1, , ]
    attr(f_lambda, "frequency") = attr(sdf.output, "frequency")[-1]
  } else {
    attr(f_lambda, "frequency") = attr(sdf.output, "frequency")
  }

  return(f_lambda)
}
