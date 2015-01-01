#' @title Computes quadratic form x' A x
#' 
#' @description
#' Computes the quadratic form \eqn{\mathbf{x}' \mathbf{A} \mathbf{x}} for an
#' \eqn{n \times n} matrix \eqn{\mathbf{A}} and an \eqn{n}-dimensional vector
#'  \eqn{\mathbf{x}}, i.e., a wrapper for \code{t(x) \%*\% A \%*\% x}. 
#'  
#'  \code{quadratic_form} works with real and complex valued matrices/vectors.
#' 
#' @param mat numeric; \eqn{n \times n} matrix (real or complex).
#' @param vec numeric; \eqn{n \times 1} vector (real or complex).
#' @return 
#' A real/complex value \eqn{\mathbf{x}' \mathbf{A} \mathbf{x}}.
#' @export
#' @keywords math univar
#' @examples
#'  set.seed(1)
#'  AA <- matrix(1:4, ncol = 2)
#'  bb <- matrix(rnorm(2))
#'  t(bb) %*% AA %*% bb
#'  quadratic_form(AA, bb)
#' 
#' 

quadratic_form <- function(mat, vec) {
  # computes the quadratic form vec' * mat * vec
  
  stopifnot(is.matrix(mat))
  dim.mat <- dim(mat)
  stopifnot(dim.mat[1] == dim.mat[2])
  
  vec <- matrix(vec)
  
  if (nrow(vec) != dim.mat[1]) {
    stop("Dimension of vector must match the matrix dimension.")
  }
  
  qp <- c(t(vec) %*% mat %*% vec)
  
  # convert to real if imaginary part is 0
  if (round(Im(qp), 6) == 0) {
    qp <- Re(qp)
  }
  names(qp) <- NULL
  return(qp)
} 