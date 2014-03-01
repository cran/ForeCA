#'@title whitens multivariate data
#'
#'@description
#'This function transforms a multivariate signal \eqn{\mathbf{X}} with mean 
#'\eqn{\boldsymbol \mu_X} and covariance matrix \eqn{\Sigma_{X}} to a \emph{whitened} 
#'signal \eqn{\mathbf{U}} with mean \eqn{\boldsymbol 0} and \eqn{\Sigma_U = I_n}.
#' Thus it makes the signal contemporaneously uncorrelated. 
#'
#'Same as projecting it on the principal components (PCs) and then scaling 
#'all PCs to unit variance.
#'
#'@details
#'The whitening is achieved by pre-multiplying \eqn{\mathbf{X}} with the inverse 
#'of the \emph{square root} of \eqn{\Sigma_{X}} using the eigen-decomposition
#'\deqn{
#'\Sigma_{X} = V' \Lambda V,
#'}
#'where \eqn{\Lambda} is an \eqn{n \times n} matrix with the eigenvalues 
#'\eqn{\lambda_1, \ldots, \lambda_n} in the diagonal.  The 
#'\emph{inverse square root} is defined as 
#'\eqn{\Sigma_{X}^{-1/2} = V' \Lambda^{-1/2}}, where 
#'\eqn{\Lambda^{-1/2} = diag(\lambda_1^{-1/2}, \ldots, \lambda_n^{-1/2})}.
#'
#'
#'@param data \eqn{T \times n} array
#'@return 
#'A list with the whitened data, the transformation, and other useful quantities
#'@export
#'@keywords manip
#'@examples
#'
#'set.seed(1)
#'XX <- matrix(rnorm(100), ncol = 2)
#'cov(XX)
#'UU <- whiten(XX)$U
#'cov(UU)

whiten <- function(data) {
  
  nseries <- ncol(data)
  
  out <- list(center = colMeans(data),
              scale = apply(data, 2, stats::sd))      
  
  if (identical(round(cor(data), 6), diag(1, nseries))) {
    # already uncorrelated
    UU <- data
    out <- c(out,
             list(Sigma0.5inv = diag(1, nseries),
                  Sigma0.5 = diag(1, nseries),
                  V = diag(1, nseries),
                  U = diag(1, nseries),
                  lambda = rep(1, nseries)))
  } else {
  
    Sigma.X <- cov(data)
    EE <- eigen(Sigma.X, symmetric = TRUE)
    VV <- EE$vectors
    UU <- sweep(data %*% VV, 2, sqrt(EE$values), FUN = "/")
    
    out <- c(out,
             list(U = UU,
                  V = VV,
                  Sigma0.5inv = VV %*% diag(1/sqrt(EE$values)),
                  Sigma0.5 = diag(sqrt(EE$values)) %*% t(VV),
                  lambda = EE$values))
  }
  return(out)
} 
