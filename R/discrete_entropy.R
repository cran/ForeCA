#'@title Shannon entropy for discrete pmf
#'
#'@description
#'Computes the Shannon entropy \eqn{\mathcal{H}(p)} of a discrete RV \eqn{X} taking
#'values in \eqn{\lbrace x_1, \ldots, x_n \rbrace} with probability
#'mass function (pmf) \eqn{Prob(X = x_i) \in \lbrace p_1, \ldots, p_n \rbrace}.
#'
#'@details
#'\code{\link{discrete_entropy}} uses a plug-in estimator (\code{method = ''MLE''}): 
#'\deqn{ 
#'\widehat{H}(p) = - \sum_{i=1}^{n} \widehat{p}_i \log \widehat{p}_i. 
#'}
#'
#'The smoothed estimator (\code{method = 'smoothed'}) adds a uniform distribution
#' \eqn{1/n} to every probability (and then normalizes again):
#'\deqn{ 
#'\widehat{p}_i \leftarrow \widehat{p_i} + \alpha \cdot \frac{1}{n}, \quad i=1, \ldots, n,
#'}
#' where \eqn{\alpha} is the smoothing weight parameter.
#'
#'@param probs probabilities (empirical frequencies). Must add up to \eqn{1}.
#'@param base see \code{base} argument of \code{\link{continuous_entropy}}
#'@param method method to estimate entropy; see Details below.
#'@param threshold frequencies below \code{threshold} are set to \eqn{0};
#'default \code{threshold = 0}.
#'@param smoothing.weight how much weight does the uniform smoothing term get?
#'default \code{1}.
#'@return 
#'Entropy (a non-negative real value).
#'@keywords math univar
#'@seealso \code{\link{continuous_entropy}}
#'@export
#'@examples
#'
#'pmfs <- rexp(5)
#'pmfs <- sort(pmfs / sum(pmfs))
#'
#'unif.distr <- rep(1/length(pmfs), length(pmfs))
#'
#'matplot(cbind(pmfs, unif.distr), pch = 19, ylab = "P(X = k)", xlab = "k")
#'matlines(cbind(pmfs, unif.distr))
#'legend("topleft", c("non-uniform","uniform"), pch = 19, lty = 1:2, col = 1:2)
#'
#'discrete_entropy(pmfs)
#' # uniform has largest entropy among all bounded distributions (here = log(5))
#'discrete_entropy(unif.distr)
#'# no uncertainty if one element occurs with probability 1
#'discrete_entropy(c(1,0,0,0,0)) 
#'
discrete_entropy <- function(probs, base = 2, method = c("MLE", "smoothed"), 
                             threshold = 0, smoothing.weight = 1) {
  
  stopifnot(all(round(probs, 6) >= 0))
  stopifnot(round(sum(probs), 6) == 1)
  
  method <- match.arg(method)
  if (threshold > 0) {
    probs[probs < threshold] <- 0
  }
  
  switch(method,
         smoothed = {
           probs <- probs + smoothing.weight * (1 / length(probs))
           probs <- probs / sum(probs)
         },
         MLE = {
           if (any(probs == 0)) {
             probs <- probs[probs != 0]
           }
         })

  entropy.eval <- -sum(probs * log(probs, base = base))
  attr(entropy.eval, "base") <- as.character(base)
  return(entropy.eval)
} 