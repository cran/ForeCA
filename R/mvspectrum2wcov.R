#'@title Compute (weighted) covariance matrix from frequency spectrum

#'@description
#'Estimates the (weighted) covariance matrix from the frequency spectrum 
#'(see Details).
#'
#'@details
#'The covariance matrix of a multivariate time series equals the average
#'of the power spectrum
#'\deqn{
#'\Sigma_{X} = \int_{-\pi}^{\pi} S_{X}(\lambda) d \lambda.
#'}
#'
#'A generalized covariance matrix can be obtained using a weighted average
#'\deqn{
#'\tilde{\Sigma}_X = \int_{-\pi}^{\pi} K(\lambda) S_{X}(\lambda) d \lambda,
#'}
#'where \eqn{K(\lambda)} is a kernel symmetric around \eqn{0} which integrates to
#'\eqn{1}. This allows one to remove or amplify specific frequencies in the 
#'covariance matrix estimate.
#'
#'For ForeCA \code{mvspectrum2wcov} is especially important as we use
#'\deqn{
#'K(\lambda) = -\log f_y(\lambda),
#'}
#'as the \emph{weights} (their average is not \eqn{1}!).
#'
#'@param mvspectrum.output output of class \code{"mvspectrum"}
#'@param weights a sequence of weights for each frequency. By default uses 
#'weights that average out to \code{1}.
#'@return 
#'An \eqn{n \times n} matrix.
#'@seealso \code{\link{mvspectrum}}
#'@keywords ts
#'@export
#'@examples
#'
#'# use SDF first and then SDF2mvspectrum
#'set.seed(1)
#'nn <- 20
#'YY <- cbind(rnorm(nn), arima.sim(n = nn, list(ar = 0.9)), rnorm(nn))
#'XX <- YY %*% matrix(rnorm(9), ncol = 3)
#'XX <- scale(XX, scale = FALSE, center = TRUE)
#'
#'# sample estimate of covariance matrix
#'Sigma.hat <- cov(XX)
#'dimnames(Sigma.hat) <- NULL
#'
#'# using the frequency spectrum
#'SS <- mvspectrum(XX, "wosa")
#'Sigma.hat.freq <- mvspectrum2wcov(SS)
#'
#'plot(c(Sigma.hat/Sigma.hat.freq))
#'abline(h = 1)
#'
#'image(Sigma.hat)
#'image(Sigma.hat.freq)
#'image(Sigma.hat / Sigma.hat.freq)
#'

mvspectrum2wcov <- function(mvspectrum.output, weights = 1) {
  nseries <- 1
  if (!is.null(dim(mvspectrum.output))) {
    nseries <- dim(mvspectrum.output)[2]
  }
  
  if (nseries == 1) {
    mvspectrum.output <- array(mvspectrum.output, 
                               c(length(mvspectrum.output), 1, 1))
  }
  
  # change the function to be integrated if weights != 1
  if (any(weights != 1)) {
    mvspectrum.output <- sweep(mvspectrum.output, 1, weights, "*")
  } 
  
  weighted.cov <- apply(mvspectrum.output, 2:3, mean)
  
  weighted.cov <- 2 * Re(weighted.cov)  # same as: weighted.cov + Conj(weighted.cov)
  return(weighted.cov)
} 