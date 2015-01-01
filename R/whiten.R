#' @title whitens multivariate data
#' 
#' @description
#' \code{whiten} transforms a multivariate K-dimensional signal \eqn{\mathbf{X}} with mean 
#' \eqn{\boldsymbol \mu_X} and covariance matrix \eqn{\Sigma_{X}} to a \emph{whitened} 
#' signal \eqn{\mathbf{U}} with mean \eqn{\boldsymbol 0} and \eqn{\Sigma_U = I_K}.
#' Thus it centers the signal and makes it contemporaneously uncorrelated. 
#' See Details.
#' 
#' @details
#' \code{whiten} uses zero component analysis (ZCA) (aka zero-phase whitening filters)
#' to whiten the data; i.e., it uses the
#' inverse square root of the covariance matrix of \eqn{\mathbf{X}} (see 
#' \code{\link{sqrt_matrix}}) as the whitening transformation.
#' This means that on top of PCA, the uncorrelated principal components are 
#' back-transformed to the original space using the
#' transpose of the eigenvectors. The advantage is that this makes them comparable
#' to the original \eqn{\mathbf{X}}.  See References for details.
#' 
#' @param data \eqn{n \times K} array representing \code{n} observations of 
#' \code{K} variables.
#' @return 
#' \code{whiten} returns a list with the whitened data, the transformation, 
#' and other useful quantities
#' @references
#' See appendix in \url{www.cs.toronto.edu/~kriz/learning-features-2009-TR.pdf}.
#' 
#' See \url{ufldl.stanford.edu/wiki/index.php/Implementing_PCA/Whitening}.
#' @export
#' @keywords manip
#' @examples
#' 
#' XX <- matrix(rnorm(100), ncol = 2) %*% matrix(runif(4), ncol = 2)
#' cov(XX)
#' UU <- whiten(XX)$U
#' cov(UU)

whiten <- function(data) {

  if (is.null(dim(data))) {
    data <- cbind(data)
  }
  num.series <- ncol(data)
  out <- list(center = colMeans(data),
              scale = apply(data, 2, stats::sd))
  data.centered <- sweep(data, 2, out$center, FUN = "-")
  
  if (num.series == 1) {
    # univariate
    out <- c(out,
             list(values = rep(var(data), num.series)))
    out <- c(out, 
             list(whitening = matrix(1 / sqrt(out$values)),
                  dewhitening = matrix(sqrt(out$values))))
    out$U <- data.centered %*% out$whitening
  } else {
    if (isTRUE(all.equal(target = diag(1, num.series), current = cov(data)))) {
      # already uncorrelated and variance 1
      out <- c(out,
               list(values = rep(1, num.series),
                    whitening = diag(1, num.series),
                    dewhitening = diag(1, num.series),
                    U = data.centered))
    } else {
      result.sqrt <- sqrt_matrix(cov(data), return.sqrt.only = FALSE, symmetric = TRUE)
      out <- c(out,
               list(values = result.sqrt$values,
                    whitening = result.sqrt$sqrt.inverse,
                    dewhitening = result.sqrt$sqrt))
      out$U <- data.centered %*% out$whitening
    }
  }
  if (num.series == 1) {
    out$U <- c(out$U)  
  }
  return(out)
} 

#' @rdname whiten
#' @export
#' @description
#' \code{check_whitened} checks if data has been whitened; i.e., if it has
#' zero mean, unit variance, and is uncorrelated.
#' @inheritParams whiten
#' @return
#' \code{check_whitened} throws an error if the input is not \code{\link{whiten}}ed.

check_whitened <- function(data) {
  
  if (is.null(dim(data))) {
    num.vars <- 1
    data <- cbind(data)
  } else {
    num.vars <- ncol(data)
  }  
  center <- colMeans(data)
  sd.data <- apply(data, 2, sd)
  
  equals <- list()
  equals[["zero-mean"]] <- all.equal(target = rep(0, num.vars), 
                                     current = center,
                                     check.names = FALSE,
                                     check.attributes = FALSE)
  equals[["unit-variance"]] <- all.equal(target = rep(1, num.vars), 
                                         current = sd.data,
                                         check.names = FALSE,
                                         check.attributes = FALSE)
  equals[["uncorrelated"]] <- all.equal(target = diag(1, num.vars),
                                        current = cor(data),
                                        check.names = FALSE,
                                        check.attributes = FALSE)
  
  errors <- c("zero-mean" = "Each variable must have mean 0.",
              "unit-variance" = "Each variable must have variance 1.",
              "uncorrelated" = "Data must be uncorrelated.")
  
  ind.errors <- !sapply(equals, isTRUE)  
  if (any(ind.errors) > 0) {
    for (ii in seq_along(ind.errors)) {
      if (ind.errors[ii]) {
        cat(names(equals)[ii], ": ", equals[[ii]], "\n", sep = "")
      }
    }
    stop("Data must be whitened:\n \t ", 
         paste(errors[ind.errors], collapse = "\n \t "))
  }
}

#' @rdname whiten
#' @description
#' \code{sqrt_matrix} computes the square root \eqn{\mathbf{B}} of a square matrix 
#' \eqn{\mathbf{A}}. The matrix \eqn{\mathbf{B}} satisfies 
#' \eqn{\mathbf{B} \mathbf{B} = \mathbf{A}}. 
#' @param mat a square \eqn{K \times K} matrix.
#' @param return.sqrt.only logical; if \code{TRUE} (default) it returns only the square root matrix;
#' if \code{FALSE} it returns other auxiliary results (eigenvectors and 
#' eigenvalues, and inverse of the square root matrix).
#' @param symmetric logical; if \code{TRUE} the \code{eigen}-solver assumes
#' that the matrix is symmetric (which makes it much faster).  This is in particular
#' useful for a covariance matrix (which is used in \code{whiten}). Default: \code{FALSE}.
#' @details
#' The \emph{square root} of a quadratic \eqn{n \times n} matrix \eqn{\mathbf{A}} 
#' can be computed by using the eigen-decomposition of \eqn{\mathbf{A}}
#' \deqn{
#'  \mathbf{A} = \mathbf{V} \Lambda \mathbf{V}',
#' }
#' where \eqn{\Lambda} is an \eqn{n \times n} matrix with the eigenvalues 
#' \eqn{\lambda_1, \ldots, \lambda_n} in the diagonal.  
#' The square root is simply \eqn{\mathbf{B} = \mathbf{V} \Lambda^{1/2} \mathbf{V}'} where 
#' \eqn{\Lambda^{1/2} = diag(\lambda_1^{1/2}, \ldots, \lambda_n^{1/2})}.
#' 
#' Similarly, the \emph{inverse square root} is defined as 
#' \eqn{\mathbf{A}^{-1/2} = \mathbf{V} \Lambda^{-1/2} \mathbf{V}'}, where 
#' \eqn{\Lambda^{-1/2} = diag(\lambda_1^{-1/2}, \ldots, \lambda_n^{-1/2})} 
#' (provided that \eqn{\lambda_i \neq 0}).
#' @return
#' \code{sqrt_matrix} returns an \eqn{n \times n}  matrix.  If \eqn{\mathbf{A}} 
#' is not semi-positive definite it returns a complex-valued \eqn{\mathbf{B}} 
#'  (since square root of negative eigenvalues are complex).
#' 
#' If \code{return.sqrt.only = FALSE} then it returns a list with:
#' 
#' \item{values}{eigenvalues of \eqn{\mathbf{A}},}
#' \item{vectors}{eigenvectors of \eqn{\mathbf{A}},}
#' \item{sqrt}{square root matrix \eqn{\mathbf{B}},}
#' \item{sqrt.inverse}{inverse of \eqn{\mathbf{B}}.}
#' @keywords math
#' @export
#' 

sqrt_matrix <- function(mat, return.sqrt.only = TRUE, symmetric = FALSE) {

  stopifnot(is.matrix(mat),
            ncol(mat) == nrow(mat))
  
  mat.eig <- eigen(mat, symmetric = symmetric)
  
  lambdas <- mat.eig$values
  if (any(lambdas < 0)) {
    lambdas[lambdas < 0] <- lambdas[lambdas < 0] + 0i
  }
  
  if (any(lambdas == 0)) {
    warning("The matrix does not have full rank; some eigenvalues are 0.")
  }
  
  if (length(lambdas) > 1) {
    diag.mat <- diag(sqrt(lambdas))
  } else {
    diag.mat <- diag(1) * sqrt(lambdas)
  }
  sqrt.mat <- as.matrix(mat.eig$vector) %*% diag.mat %*% t(mat.eig$vector)
  if (return.sqrt.only) {
    return(sqrt.mat)
  } else {
    diag(diag.mat) <- 1 / diag(diag.mat)
    sqrt.mat.inverse <- mat.eig$vector %*% diag.mat %*% t(mat.eig$vector)
    out <- list(values = lambdas,
                vectors = mat.eig$vector,
                sqrt = sqrt.mat,
                sqrt.inverse = sqrt.mat.inverse)
    return(out)
  }
}

