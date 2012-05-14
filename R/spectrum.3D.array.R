spectrum.3D.array <-
function(series, method = "wosa", frequency_0 = TRUE, ...){
  sdf = SDF(series, method = method, ...)
  return(SDF2spectrum.3D.array(sdf.output = sdf, frequency_0 = frequency_0))
}
