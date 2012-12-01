#'@title Shannon entropy for a continuous pdf
#'
#'@description
#'Computes the Shannon entropy \eqn{H(p)} for a continuous probability density function
#'(pdf) \eqn{p(x)}.
#'
#'@details
#'The Shannon entropy of a continuous RV \eqn{X \sim p(x)} is defined as
#'\deqn{
#'H(p) = -\int_{-\infty}^{\infty} p(x) \log p(x) d x
#'}
#'
#'Contrary to the entropy of a discrete RV, continuous RVs can
#'have negative entropy (see Examples).
#'
#'@param pdf R function (\code{function(x) ...}) for the pdf \eqn{p(x)} 
#'of a RV \eqn{X \sim p(x)}
#'@param base logarithm base; default \code{base = 2}. Entropy is measured in
#'``nats'' for \code{base = e}; in ``bits'' if \code{base = 2}.
#'@param a lower integration limit
#'@param b upper integration limit
#'@return 
#'Entropy (real). 
#'@keywords math univar
#'@seealso \code{\link{discrete_entropy}}
#'@export
#'@examples
#'
#'my_density = function(x){
#'  dunif(x, -1, 1)
#'}
#'
#'\dontrun{
#'# an error since limits of integration were not specified
#'continuous_entropy(my_density)
#'}
#'continuous_entropy(my_density, -1, 1) # = log(b - a)
#'
#'# entropy of the U(a,b) = log(b-a). thus not necessarily positive anymore, e.g.
#'continuous_entropy(function(x) dunif(x, 0, 0.5), 0, 0.5) # log2(0.5)
#'
#'
continuous_entropy = function(pdf = function(x) return(dunif(x, 0, 1)), 
                              a = NULL, b = NULL, base = 2){
  aux = function(xx){
    - pdf(xx) * log(pdf(xx), base = base)
  }
  
  if (is.null(a)){
    stop("You must provide the lower limit for which p(x)>0.")
  }
  if (is.null(b)){
    stop("You must provide the upper limit for which p(x)>0.")
  }
  
  HH = integrate(aux, a, b)$value
  attr(HH, "base") = as.character(base)
  return(HH)
}




