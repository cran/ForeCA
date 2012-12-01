#' @title ForeCA EM auxiliary functions
#' @name ForeCA.EM-aux
#' @aliases ForeCA.EM.E_step ForeCA.EM.M_step ForeCA.EM.h
#' @description 
#' The EM-like algorithm needs three auxiliary functions:
#' 
NULL

#' @rdname ForeCA.EM-aux
#' @description
#' \code{ForeCA.EM.E_step} computes the spectral density given a weight-vector
#' \eqn{\mathbf{w}_i}.
#' @keywords manip
#' @param f_U an object of class \code{"mvspectrum"} with 
#' \code{normalized = TRUE}
#' @param weights a weight vector for \eqn{y_t=\mathbf{w}' \mathbf{X}_t}. 
#' Must have unit norm in \eqn{L_2}. 
#' @export
#' @return
#' \code{ForeCA.EM.E_step} returns a vector containing the normalized spectral 
#' density (normalized such that its \code{mean} equals \eqn{0.5} - since 
#' \code{f_U} only contains positive frequencies and is symmetric).
#' @examples
#' 
#' XX = diff(log(EuStockMarkets)) * 100
#' ff = mvspectrum(XX, 'wosa', normalize = TRUE)
#' 
#' ww0 = matrix(rnorm(ncol(XX)))
#' ww0 = ww0 / norm(ww0, 'F')
#' 
#' f_ww0 = ForeCA.EM.E_step(ff, ww0)
#' 

ForeCA.EM.E_step <- function(f_U, weights = NULL) {
  weights <- weights/base::norm(as.matrix(weights), "F")
  ff <- Re(apply(f_U, 1, quadratic_form, vec = weights))
  ff[ff < 0] <- 1/length(ff)^2
  
  invisible(normalize_mvspectrum(ff))
} 


#' @rdname ForeCA.EM-aux
#' @description
#' \code{ForeCA.EM.M_step} computes the minimizing eigenvector 
#' (\eqn{\rightarrow \widehat{\mathbf{w}}_{i+1}}) of the weighted
#' covariance matrix, where the weights equal the negative logarithm of the 
#' spectral density at the current \eqn{\widehat{\mathbf{w}}_i}.
#' @keywords manip
#' @param base logarithm base; if \code{NULL} it sets it automatically to \eqn{T}, 
#' such that the resulting entropy estimate is bounded above by \eqn{1} (it is always
#' bounded below by \eqn{0}).
#' @param minimize indicator; if \code{TRUE} it returns the eigenvector
#' correpsonding to the smallest eigenvalue; otherwise to the largest eigenvalue.
#' @param Sigma_X optional; covariance matrix of \eqn{\mathbf{X}}
#' @export
#' @return
#' \code{ForeCA.EM.M_step} returns a list with three elements:
#' \describe{
#'    \item{\code{matrix}:}{is the weighted covariance matrix 
#'          (positive semi-definite), where the weights are the negative 
#'          log of the spectral density}
#'    \item{\code{vector}:}{minimizing (or maximizing if 
#'          \code{minimize = FALSE}) eigenvector of that matrix}
#'    \item{\code{value}:}{corresponding eigenvalue}
#' }
#' @examples
#' 
#' onestep = ForeCA.EM.M_step(ff, f_ww0)
#' image(onestep$matrix)
#' \dontrun{
#' # if you have 'LICORS' package installed, better use this:
#' image2(onestep$matrix)
#' }
#' ww1 = onestep$vector
#' f_ww1 = ForeCA.EM.E_step(ff, ww1)
#' 
#' par(mar = c(4, 4, 1, 1), mfcol= c(1,2))
#' matplot(seq(0, pi, length = length(f_ww0)), cbind(f_ww0, f_ww1), type = "l", lwd =2, xlab = "frequencies")
#' plot(f_ww0, f_ww1, pch = ".", cex = 3, xlab = "iteration 0", ylab = "iteration 1")
#' 
#' 

ForeCA.EM.M_step <- function(f_U, f_current = NULL, base = NULL, 
                             minimize = TRUE, 
                             Sigma_X = NULL) {
  
  TT <- length(f_current)
  nseries <- ncol(f_U)
  #f_current <- f_current #/ sum(f_current) / 2
  
  if (round(mean(f_current), 6) != 0.5){
    stop("The spectral density in 'f_current' is not correctly normalized. Its
         mean must equal 0.5 (since it's symmetric around 0 we only have 
         the spectral density for omega > 0.")
  }
  
  if (!identical(round(mvspectrum2wcov(f_U), 6), diag(1, nseries))){
    stop("The multivariate spectrum estimate 'f_U' must be properly normalized.
         It must average (mean) to a diagonal matrix with 0.5 in the diagonal.")
  }
  if (is.null(base)){
    base <- 2 * TT
  }
  SS <- mvspectrum2wcov(f_U / TT, -log(f_current / TT, base = base))
  
  if (!is.null(Sigma_X)){
    SS <- solve(Sigma_X %*% SS)
  }
  
  EE <- eigen(SS, symmetric = TRUE)
  if (minimize) {
    sel <- which.min(EE$values)
  } else {
    sel <- which.max(EE$values)
  }
  weights <- EE$vector[, sel]
  value <- EE$values[sel]
  
  weights <- weights/sign(weights[which.max(abs(weights))])
  weights <- weights/base::norm(as.matrix(weights), "F")
  
  out <- list()
  out$matrix <- SS
  out$vector <- weights
  out$value <- value
  return(out)
} 

#' @rdname ForeCA.EM-aux
#' @description
#' \code{ForeCA.EM.h} evaluates the entropy of the spectral density as a function
#' of \eqn{\mathbf{w}}. This is the objective funcion that should be minimized.
#' 
#' @keywords manip
#' @param weights_new weights \eqn{\widehat{\mathbf{w}}_{i+1}} of the new iteration (i+1)
#' @param f_current spectral density of the current estimate 
#' \eqn{\widehat{\mathbf{w}}_i} (required for \code{ForeCA.EM.M_step}; 
#' optional for \code{ForeCA.EM.h}).
#' @param weights_current weights \eqn{\widehat{\mathbf{w}}_{i}} of the current iteration (i)
#' @return 
#' \code{ForeCA.EM.h} returns (see References for details):
#' \itemize{
#'    \item the negative entropy (non-negative real) if 
#'          \code{weights_new = weights_current}
#'    \item an upper bound of that entropy for \code{weights_new} if 
#'          \code{weights_new != weights_current}
#' }
#' @export
#' @examples
#' 
#' ForeCA.EM.h(ww0, ff)       # iteration 0
#' ForeCA.EM.h(ww1, ff, ww0)  # min eigenvalue inequality
#' ForeCA.EM.h(ww1, ff)       # KL divergence inequality
#' onestep$value
#' 
#' Omega(spectrum_estimate = f_ww0) / 100 + ForeCA.EM.h(ww0, ff)
#' Omega(spectrum_estimate = f_ww1) / 100 + ForeCA.EM.h(ww1, ff)
#' 

ForeCA.EM.h <- function(weights_new, f_U, 
                        weights_current = weights_new, f_current = NULL, 
                        base = NULL) {
  # short as quadratic_form(apply(f_U*log(weights_current), 2:3, sum),
  # weights_new)
  if (is.null(f_current)) {
    f_current <- ForeCA.EM.E_step(f_U, weights_current)  #Re(apply(f_U, 1, quadratic_form, weights_current))
  }
  if (round(mean(f_current), 6) != 0.5){
    stop("The spectral density in 'f_current' is not correctly normalized. Its
         mean must equal 0.5 (since it's symmetric around 0 we only have 
         the spectral density for omega > 0.")
  }
  
  nseries <- length(weights_new)
  
  if (!identical(round(mvspectrum2wcov(f_U), 6), diag(1, nseries))){
    stop("The multivariate spectrum estimate 'f_U' must be properly normalized.
         It must average (mean) to a diagonal matrix with 0.5 in the diagonal.")
  }
  
  TT <- length(f_current)
  
  if (is.null(base)) {
    base <- 2 * length(f_current)
  }
  
  wCovMatrix <- mvspectrum2wcov(f_U/TT, weights = -log(f_current / TT, base = base))
  
  return(quadratic_form(wCovMatrix, weights_new))
} 


#' @rdname ForeCA.EM-aux
#' @description
#' \code{ForeCA.EM.init_weightvector} returns an initialization vector 
#' \eqn{\mathbf{w}_0 \in R^k}, with \eqn{L_2} norm \eqn{1}, for
#' the EM-like algorithm in \code{\link{ForeCA.EM.one_weightvector}}.  
#' They are all based on quickly computable heuristics
#' than approximate a forecastable signal. See \code{method} argument and Details
#' below. 
#' 
#' @keywords manip
#' @param UU uncorrelated multivariate time series
#' @param method string indicating the method to estimate the initial 
#' weightvector; default \code{"SFA"}. See Details.
#' @param seed seed to use for random initialization; if \code{NULL}, it sets
#' a random seed
#' @param lag integer; lag for the difference operator; default \code{lag=1}.
#' @details
#' The \code{method} argument determines what heuristics should be used to 
#' get a good starting vector \eqn{\mathbf{w}_0}. This vector is length \eqn{k}
#' and has unit norm in \eqn{L_2}: 
#' 
#' \describe{
#'  \item{max}{uses the vector with all \eqn{0}s, but a \eqn{1} at the position
#'  of the maximum forecastability of the series in \code{UU}.}
#'  \item{SFA_slow}{uses the first eigenvector of SFA (slowest signal).}
#'  \item{SFA_fast}{uses the last eigenvector of SFA (fastes signal).}
#'  \item{SFA}{checks both first and last, and chooses the one with higher
#'             forecastability.}
#'  \item{cauchy}{random starts using \code{rcauchy(k)}}
#'  \item{unif}{random starts using \code{runif(k, -1, 1)}}
#'  \item{norm}{random starts using \code{runif(k, 0, 1)}}
#' }
#' 
#' Slow Feature Analysis (SFA; see References) finds \emph{slow} signals. For 
#' details see the References; for us it's important that SFA is equivalent to
#' finding the signal with largest lag \eqn{1} autocorrelation.  This is not
#' necessarily the most forecastable, but a good start.
#' 
#' The disadvantage of SFA for forecasting is that e.g. white noise (WN) 
#' is ranked higher than an AR(1) with negative autocorrelation.  While it
#' is true that WN is slower, it is not more forecastable.  Thus we are also 
#' interested in the fastest signal, i.e. the last eigenvector.
#'                  determine a good starting vector. The slowest feature
#'                  corresponds to maximizing the lag \eqn{1} correlation
#' @return 
#' \code{ForeCA.EM.init_weightvector} returns a good starting vector 
#' \eqn{\mathbf{w}_0} for \code{\link{ForeCA.EM.one_weightvector}}.
#' @references
#' Laurenz Wiskott and Terrence J. Sejnowski (2002). 
#' \dQuote{Slow Feature Analysis: Unsupervised Learning of Invariances}, 
#' Neural Computation 14:4, 715-770.

ForeCA.EM.init_weightvector <- function(UU, f_U, method="SFA", seed=NULL, 
                                        lag = 1) {
  kk <- ncol(UU)
  
  if (!identical(round(cov(UU), 3), diag(1, kk))){
    stop("Time series must be uncorrelated and have unit variance.")
  }
  
  if (is.null(seed)){
    seed <- sample.int(1e6, 1)
  }
  if (method == "cauchy"){
    ww0 <- rcauchy(kk)
  } else if (method == 'unif'){
    ww0 <- runif(kk, -1, 1)
  } else if (method == "norm"){
    ww0 <- rnorm(kk)
  } else if (method == "max"){
    ww0 <- rep(0, kk)
    omega_temp <- apply(2*Re(apply(f_U, 1, diag)), 1, spectral_entropy, 
                        threshold = 0)
    ww0[which.max(omega_temp)] <- 1
  } else if (any(method == c("SFA", "SFA_slow", "SFA_fast"))){
    Sigma_Delta_U <- cov(diff(UU, lag = lag))
    EE <- eigen(Sigma_Delta_U)
    ww_slow <- EE$vectors[, kk]
    ww_fast <- EE$vectors[, 1]
    
    ww0 <- ww_slow
    if (method == "SFA_slow"){
      # done
    } else if (method == "SFA_fast"){
      ww0 <- ww_fast
    } else if (method == "SFA"){
      omega_slow <- Omega(spectrum_estimate = ForeCA.EM.E_step(f_U, ww_slow))
      omega_fast <- Omega(spectrum_estimate = ForeCA.EM.E_step(f_U, ww_fast))
      
      if (omega_fast > omega_slow){
        ww0 <- ww_fast
      }
    }
    
  }
  ww0 <- ww0 / base::norm(ww0, "2")
  
  return(ww0)
} 


