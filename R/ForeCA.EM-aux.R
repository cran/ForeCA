#' @title ForeCA EM auxiliary functions
#' @name foreca.EM-aux
#' @aliases foreca.EM.E_step foreca.EM.M_step foreca.EM.h
#' @description 
#' The EM-like algorithm relies on several auxiliary functions:
#' 
NULL

#' @rdname foreca.EM-aux
#' @description
#' \code{foreca.EM.E_step} computes the spectral density of 
#' \eqn{y_t=\mathbf{w}' \mathbf{X}_t} given the weightvector \eqn{\mathbf{w}}.
#' @keywords manip
#' @param f.U an object of class \code{"mvspectrum"} with 
#' \code{normalized = TRUE}
#' @param weightvector a weight vector for \eqn{y_t=\mathbf{w}' \mathbf{X}_t}. 
#' Must have unit norm in \eqn{L^2}. 
#' @export
#' @return
#' \code{foreca.EM.E_step} returns a vector containing the normalized spectral 
#' density (normalized such that its \code{mean} equals \eqn{0.5} - since 
#' \code{f.U} only contains positive frequencies and is symmetric).
#' @examples
#' 
#' XX <- diff(log(EuStockMarkets)) * 100
#' ff <- mvspectrum(XX, 'wosa', normalize = TRUE)
#' 
#' ww0 <- matrix(rnorm(ncol(XX)))
#' ww0 <- ww0 / norm(ww0, 'F')
#' 
#' f.ww0 <- foreca.EM.E_step(ff, ww0)
#' 

foreca.EM.E_step <- function(f.U, weightvector) {
  
  stopifnot(round(base::norm(as.matrix(weightvector), "F"), 6) == 1)
  
  spec.dens.est <- Re(apply(f.U, 1, quadratic_form, vec = weightvector))
  spec.dens.est[spec.dens.est < 0] <- 1 / length(spec.dens.est)^2
  
  return(normalize_mvspectrum(spec.dens.est))
} 


#' @rdname foreca.EM-aux
#' @description
#' \code{foreca.EM.M_step} computes the minimizing eigenvector 
#' (\eqn{\rightarrow \widehat{\mathbf{w}}_{i+1}}) of the weighted
#' covariance matrix, where the weights equal the negative logarithm of the 
#' spectral density at the current \eqn{\widehat{\mathbf{w}}_i}.
#' @keywords manip
#' @param base logarithm base; if \code{NULL} it sets it automatically to \eqn{T}, 
#' such that the resulting discrete entropy estimate is bounded above by \eqn{1} 
#' (it is always bounded below by \eqn{0}).
#' @param minimize logical; if \code{TRUE} it returns the eigenvector
#' corresponding to the smallest eigenvalue; otherwise to the largest eigenvalue.
#' @param Sigma.X optional; covariance matrix of \eqn{\mathbf{X}}
#' @export
#' @return
#' \code{foreca.EM.M_step} returns a list with three elements:
#' \describe{
#'    \item{\code{matrix}:}{is the weighted covariance matrix 
#'          (positive semi-definite), where the weights are the negative 
#'          log of the spectral density,}
#'    \item{\code{vector}:}{minimizing (or maximizing if 
#'          \code{minimize = FALSE}) eigenvector of \code{matrix},}
#'    \item{\code{value}:}{corresponding eigenvalue.}
#' }
#' @examples
#' 
#' onestep <- foreca.EM.M_step(ff, f.ww0)
#' image(onestep$matrix)
#' \dontrun{
#' require(LICORS)
#' # if you have 'LICORS' package installed, better use this:
#' image2(onestep$matrix)
#' }
#' ww1 <- onestep$vector
#' f.ww1 <- foreca.EM.E_step(ff, ww1)
#' 
#' layout(matrix(1:2, ncol = 2))
#' #par(mar = c(4, 4, 1, 1), mfcol= c(1,2))
#' matplot(seq(0, pi, length = length(f.ww0)), cbind(f.ww0, f.ww1), 
#'         type = "l", lwd =2, xlab = "omega_j", ylab = "f(omega_j)")
#' plot(f.ww0, f.ww1, pch = ".", cex = 3, xlab = "iteration 0", 
#'      ylab = "iteration 1", main = "Spectral density")
#' abline(0, 1, col = 'blue', lty = 2)
#' 
#' 

foreca.EM.M_step <- function(f.U, f.current, base = NULL, minimize = TRUE, 
                             Sigma.X = NULL) {
  
  TT <- length(f.current)
  num.series <- ncol(f.U)
  
  if (round(mean(f.current), 6) != 0.5) {
    stop(paste("The spectral density in 'f.current' is not correctly normalized.",
               "\n Its mean must equal 0.5 (since it's symmetric around 0 only",
               "the spectral density for omega > 0 is used)."))
  }
  
  if (!identical(round(mvspectrum2wcov(f.U), 6), diag(1, num.series))) {
    stop(paste("The multivariate spectrum estimate 'f.U' must be properly",
               " normalized.\n It must average (mean) to a diagonal matrix ",
               " with 0.5 in the diagonal."))
  }
  if (is.null(base)) {
    base <- 2 * TT
  }
  SS <- mvspectrum2wcov(f.U/TT, -log(f.current/TT, base = base))
  
  if (!is.null(Sigma.X)) {
    SS <- solve(Sigma.X %*% SS)
  }
  
  EE <- eigen(SS, symmetric = TRUE)
  sel <- ifelse(minimize, 
                which.min(EE$values), max(EE$values))
  
  weightvector <- EE$vector[, sel]
  value <- EE$values[sel]
  
  weightvector <- weightvector / sign(weightvector[which.max(abs(weightvector))])
  weightvector <- weightvector / base::norm(as.matrix(weightvector), "F")
  
  out <- list(matrix = SS, vector = weightvector, value = value)
  return(out)
} 

#' @rdname foreca.EM-aux
#' @description
#' \code{foreca.EM.h} evaluates the entropy of the spectral density as a function
#' of \eqn{\mathbf{w}}. This is the objective funcion that should be minimized.
#' 
#' @keywords manip
#' @param weightvector.new weightvector \eqn{\widehat{\mathbf{w}}_{i+1}} of the new 
#' iteration (i+1)
#' @param f.current spectral density of \eqn{y_t=\mathbf{w}' \mathbf{X}_t} for 
#' the current estimate \eqn{\widehat{\mathbf{w}}_i} (required for 
#' \code{foreca.EM.M_step}; optional for \code{foreca.EM.h}).
#' @param weightvector.current weightvector \eqn{\widehat{\mathbf{w}}_{i}} of the 
#' current iteration (i)
#' @return 
#' \code{foreca.EM.h} returns (see References for details):
#' \itemize{
#'    \item the negative entropy (non-negative real) if 
#'          \code{weightvector.new = weightvector.current}
#'    \item an upper bound of that entropy for \code{weightvector.new} if 
#'          \code{weightvector.new != weightvector.current}
#' }
#' @export
#' @examples
#' 
#' foreca.EM.h(ww0, ff)       # iteration 0
#' foreca.EM.h(ww1, ff, ww0)  # min eigenvalue inequality
#' foreca.EM.h(ww1, ff)       # KL divergence inequality
#' onestep$value
#' 
#' Omega(spectrum.estimate = f.ww0) / 100 + foreca.EM.h(ww0, ff)
#' Omega(spectrum.estimate = f.ww1) / 100 + foreca.EM.h(ww1, ff)
#' 

foreca.EM.h <- function(weightvector.new, f.U, 
                        weightvector.current = weightvector.new, 
                        f.current = NULL, base = NULL) {
  # short as quadratic_form(apply(f.U*log(weightvector.current), 2:3, sum),
  # weightvector.new)
  if (is.null(f.current)) {
    f.current <- foreca.EM.E_step(f.U, weightvector.current)
  }
  if (round(mean(f.current), 6) != 0.5) {
    stop(paste("The spectral density in 'f.current' is not correctly normalized.",
               "\n Its mean must equal 0.5 (since it's symmetric around 0 only",
               "spectral density for omega > 0 is required."))
  }
  
  num.series <- length(weightvector.new)
  
  if (!identical(round(mvspectrum2wcov(f.U), 6), diag(1, num.series))) {
    stop(paste("The multivariate spectrum estimate 'f.U' must be properly ",
               "normalized.\n It must average (mean) to a diagonal matrix ",
               "with 0.5 in the diagonal."))
  }
  
  TT <- length(f.current)
  
  if (is.null(base)) {
    base <- 2 * length(f.current)
  }
  
  w.cov.matrix <- 
    mvspectrum2wcov(f.U/TT, weights = -log(f.current/TT, base = base))
  
  return(quadratic_form(w.cov.matrix, weightvector.new))
} 

