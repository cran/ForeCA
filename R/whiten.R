#'@title whitens multivariate data
#'
#'@description
#'This function transforms a multivariate signal \eqn{\mathbf{X}} with mean 
#'\eqn{\mu_X} and covariance matrix \eqn{\Sigma_X} to a \emph{whitened} 
#'signal \eqn{\mathbf{U}} with mean \eqn{0} and \eqn{\Sigma_U = I_n}. That is 
#'it makes the signal contemporaneously uncorrelated. 
#'
#'Same as projecting it on the principal components (PCs) and then scaling 
#'all PCs to unit variance.
#'
#'@details
#'The whitening is achieved by pre-multiplying \eqn{\mathbf{X}} with the inverse 
#'of the \emph{square root} of \eqn{\Sigma_X} using the eigen-decomposition
#'\deqn{
#'\Sigma_X = V' \Lambda V,
#'}
#'where \eqn{\Lambda} is an \eqn{n \times n} matrix with the eigenvalues 
#'\eqn{\lambda_1, \ldots, \lambda_n} in the diagonal.  The 
#'\emph{inverse square root} is defined as 
#'\eqn{\Sigma_X^{-1/2} = V' \Lambda^{-1/2}}, where 
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
#'XX = matrix(rnorm(100), ncol = 2)
#'cov(XX)
#'UU = whiten(XX)$U
#'cov(UU)

whiten <- function(data) {
  # Compute an uncorrelated, unit-variance version of X, using the SVD form
  # of cov(XX) = Sigma_X = V' Lambda V and thus get U = V' Lambda^(-1/2) X
  
  nseries <- ncol(data)
  
  out <- list()
  
  out$center <- colMeans(data)
  out$scale <- apply(data, 2, stats::sd)
  
  if (identical(round(cor(data),6), diag(1, nseries))){
    # already uncorrelated
    UU <- data
    out$Sigma0.5inv <- diag(1, nseries)
    out$Sigma0.5 <- diag(1, nseries)
    out$V <- diag(1, nseries)
    out$U <- diag(1, nseries)
    out$lambda <- rep(1, nseries)
  } else {
    #data <- scale(data, scale = FALSE, center = TRUE)
  
    Sigma_X <- cov(data)
    EE <- eigen(Sigma_X)
    VV <- EE$vectors
    PP <- data %*% VV
    UU <- sweep(PP, 2, sqrt(EE$values), FUN = "/")
    
    out$U <- UU
    out$V <- VV
    out$Sigma0.5inv <- VV %*% diag(1/sqrt(EE$values))
    out$Sigma0.5 <- diag(sqrt(EE$values)) %*% t(VV)
    out$lambda <- EE$values
  }
  invisible(out)
} 