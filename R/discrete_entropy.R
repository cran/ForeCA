#'@title Shannon entropy for discrete pmf
#'
#'@description
#'Computes the Shannon entropy \eqn{\mathcal{H}(p)} of a discrete probability
#'mass function (pmf) \eqn{\lbrace p_1, \ldots, p_n \rbrace}.
#'
#'@details
#'\code{\link{discrete_entropy}} uses a plug-in estimator (\code{method =
#'''MLE''}): 
#'\deqn{ 
#'\widehat{H}(p) = - \sum_{i=1}^{n} \widehat{p}_i \log \widehat{p}_i. 
#'}
#'
#'The smoothed estimator (\code{method = 'smoothed'}) adds \eqn{1/n} to every 
#'probability (and then normalizes again):
#'\deqn{ 
#'\widehat{p}_i \leftarrow \widehat{p_i} + \frac{1}{n}, \quad i=1, \ldots, n. 
#'}
#'
#'@param probs probabilities (empirical frequencies)
#'@param base logarithm base; default \code{base = 2}. Entropy is measured in
#'``nats'' for \code{base = e}; in ``bits'' if \code{base = 2}.
#'@param method method to estimate entropy; see Details below.
#'@param threshold frequencies below \code{threshold} are set to \eqn{0};
#'default \code{threshold = 0}.
#'@return 
#'Entropy (a non-negative real value).
#'@keywords math univar
#'@seealso \code{\link{continuous_entropy}}
#'@export
#'@examples
#'
#'pmfs = rexp(5)
#'pmfs = sort(pmfs / sum(pmfs))
#'
#'unif.distr = rep(1/length(pmfs), length(pmfs))
#'
#'matplot(cbind(pmfs, unif.distr), pch = 19, ylab = "P(X = k)", xlab = "k")
#'matlines(cbind(pmfs, unif.distr))
#'legend("topleft", c("non-uniform","uniform"), pch = 19, lty = 1:2, col = 1:2)
#'
#'discrete_entropy(pmfs)
#'discrete_entropy(unif.distr) # uniform has largest entropy among all bounded distributions
#'discrete_entropy(c(1,0,0,0,0)) # no uncertainty if one element happens with probability 1
#'
discrete_entropy <-
function(probs, base = 2, method = "MLE", threshold = 0){
  
  probs = probs / sum(probs)
  
  if (threshold > 0) {
    probs[probs < threshold] = 0
  }
  if (method == "smoothed"){
    probs = probs + 1/length(probs)
    probs = probs / sum(probs)
    HH = - sum(probs * log(probs, base = base))
  }
  if (method == "MLE"){
    if (any(probs == 0)) {
      probs = probs[probs != 0]
    }
    HH = -sum(probs * log(probs, base = base))
  }

  attr(HH, "base") = as.character(base)
  return(HH)
}
