Omega <-
function(series, spectrum_method = "multitaper", entropy_method = "MLE", 
                 threshold = NULL, ...){
  
  series = as.matrix(series)
  # scale to zero-mean and unit-variance
  series = scale(series)
  TT = nrow(series)
  
  if (is.null(threshold)) {
    threshold = 1/(TT^(1.5))
  }
  
  if (ncol(series) > 1) {
    H_s = apply(series, 2, spectral.entropy, 
                spectrum_method = spectrum_method, base = 2, 
                entropy_method = entropy_method, 
                threshold = threshold, ...)
  } else {
    H_s = spectral.entropy(series, spectrum_method = spectrum_method, base = 2, 
                           entropy_method = entropy_method, 
                           threshold = threshold, ... )
  }
  if (entropy_method == "MLE") {
    OMEGA = 1 - H_s/(log(2*TT, base = 2))
  }
  return(c(OMEGA*100))
}
