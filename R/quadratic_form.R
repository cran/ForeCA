#'@title Computes quadratic form x' A x
#'
#'@description
#'Computes the quadratic form \eqn{\mathbf{x}' \mathbf{A} \mathbf{x}} for a 
#'\eqn{n \times n} matrix 
#'\eqn{\mathbf{A}} and an \eqn{n}-dimensional vector \eqn{\mathbf{x}}.
#'
#'Works with real and complex valued matrices/vectors.
#'
#'A wrapper for \code{t(x) \%*\% A \%*\% x}. 
#'
#'@param mat \eqn{n \times n} matrix (real or complex)
#'@param vec \eqn{n \times 1} vector (real or complex)
#'@return 
#'A real/complex value (\eqn{\mathbf{x}' \mathbf{A} \mathbf{x}}).
#'@export
#'@keywords math univar
#'@examples
#'set.seed(1)
#'AA = matrix(1:4, ncol = 2)
#'bb = matrix(rnorm(2))
#'t(bb) %*% AA %*% bb
#'quadratic_form(AA, bb)
#'
#'


quadratic_form <- function(mat, vec) {
  # computes the quadratic form vec' * mat * vec
  kk <- dim(mat)[1]
  vec <- matrix(vec)
  
  if (nrow(vec) != kk) {
    stop("Dimension of vector must match the matrix dimension.")
  }
  
  qp <- c(t(vec) %*% mat %*% vec)
  
  if (round(Im(qp), 6) == 0) {
    qp <- Re(qp)
  }
  
  return(qp)
} 