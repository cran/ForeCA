#' @title EM-like Algorithm to Estimate Optimal Transformations via ForeCA
#' @name foreca.EM
#' @aliases foreca.EM.opt_comp
#' @description 
#' 
#' \code{foreca.EM} estimates the optimal transformations to obtain forecastable
#' signals from multivariate time series \code{series} 
#' 
#' @inheritParams foreca.EM.opt_comp
#' @param n.comp number of components to be extracted
#' @param plot indicator; if \code{TRUE} a plot of the current optimal 
#' solution \eqn{\mathbf{w}_i^*} will be shown and the plot is updated for each 
#' new optimal weightvector.
#' @return
#' A list with similar output as \code{\link[stats]{princomp}}. Signals are
#' ordered from most to least forecastable. 
#' @export
#' @examples
#' \dontrun{
#' XX = diff(log(EuStockMarkets[-c(1:50),])) * 100
#' plot(ts(XX))
#' foreca = foreca.EM(XX[,1:3], n.comp = 3)
#' 
#' summary(foreca)
#' plot(foreca)
#' }

foreca.EM <- function(series, 
                      spectrum_method = "multitaper", 
                      n.comp = 2, 
                      tol = 1e-4, 
                      max_iter = 200, 
                      kernel = NULL, 
                      nstart = 10, 
                      plot = TRUE,
                      threshold = 0,
                      smoothing = FALSE, ...) {
  
  if (!is.ts(series)) {
    series <- ts(series)
  }
  
  TT <- nrow(series)
  
  PW <- whiten(series)
  UU <- PW$U
    
  kk <- ncol(UU)
  null_space <- diag(1, kk)
  
  ff_UU <- mvspectrum(UU, method = spectrum_method, normalize = TRUE, ...)
  if (!is.null(kernel)) {
    ff_UU <- sweep(ff_UU, 1, kernel(attr(ff_UU, "frequency")), FUN = "*")
    ff_UU <- normalize_mvspectrum(ff_UU)
  }
  
  out <- list()
  out$loadings <- matrix(NA, ncol = n.comp, nrow = kk)
  out$h <- rep(NA, n.comp)
  out$h.tries <- matrix(NA, ncol = n.comp, nrow = nstart)

  out$weights <- list()
  out$Omega2 <- rep(NA, n.comp)
  out$h.orig <- rep(NA, n.comp)
  out$best.tries <- rep(NA, n.comp)
  for (ii in 1:n.comp) {
    UU_ortho <- as.matrix(UU) %*% null_space
    ff_UU_ortho <- NULL
    #ff_UU_ortho <- array( t(apply(ff_UU, 1, 
    #                                   function(mat, x) t(x) %*% mat %*% x, 
    #                                                    x = null_space)), 
    #                dim = c(dim(ff_UU)[1], dim(null_space)[2], dim(null_space)[2] ) )
    #print(UU_ortho[1:2, ])
    if (ii == kk){
      weight_ortho <- 1
      out$h[ii] <- spectral_entropy(UU_ortho, 
                                    spectrum_method = spectrum_method,
                                    threshold = threshold)
      
    } else {
      foreca_ortho <- foreca.EM.opt_comp(UU_ortho, 
                                                 f_U = ff_UU_ortho,
                                                 spectrum_method = spectrum_method, 
                                                 kernel = kernel, 
                                                 nstart = nstart, 
                                                 threshold = threshold,
                                                 smoothing = smoothing)
      out$best.tries[ii] <- foreca_ortho$best.try
      out$h.orig[ii] <- foreca_ortho$h.trace[1]
      weight_ortho <- foreca_ortho$loadings
      out$weights[[ii]] <- foreca_ortho$weights
      out$Omega2[ii] <- foreca_ortho$Omega
      
      out$h[ii] <- foreca_ortho$h
      out$h.tries[, ii] <- foreca_ortho$h.tries
      
      # plot results
      if (plot) {
        plot(foreca_ortho, main = paste("Component", ii))
        Sys.sleep(0.5)
      }
    }

    out$loadings[, ii] <- null_space %*% weight_ortho
    out$loadings[, ii] <- out$loadings[, ii] / base::norm(as.matrix(out$loadings[, ii]), 
                                                        "F")
    null_space <- Null(out$loadings[, 1:ii])
  }
  
  if (any(is.na(out$loadings))) {
    out$loadings[, kk] <- Null(out$loadings[, -kk])
    out$loadings[, kk] <- out$loadings[, kk]/base::norm(as.matrix(out$loadings[, kk]), 
                                                        "F")
  }
  
  out$x <- series
  out$tol <- tol
  
  out$scores2 <- ts(as.matrix(UU) %*% out$loadings, 
                    start = start(series), end = end(series))
  
  out$loadings <- PW$Sigma0.5inv %*% out$loadings
  out$loadings <- sweep(out$loadings, 2,
                        apply(out$loadings, 2, norm, "2"), "/")
  out$Sigma0.5 <- PW$Sigma0.5
  out$Sigma0.5inv <- PW$Sigma0.5inv
  # make the sign of the loadings consistent for repeated runs
  out$loadings <- sweep(out$loadings, 2, sign(out$loadings[1, ]), "*")
  
  out$n.obs <- nrow(series)
  out$scores <- ts(scale(series, scale = FALSE, center = TRUE) %*% out$loadings, 
                   start = start(series), end = end(series))
  
  omega.scores <- Omega(out$scores, spectrum_method = spectrum_method, 
                        threshold = threshold)
  
  omega_order <- order(omega.scores, decreasing = TRUE)
  
  out$scores <- out$scores[, omega_order]
  omega.scores <- omega.scores[omega_order]
  out$loadings <- out$loadings[, omega_order]
  
  colnames(out$scores) <- paste("ForeC", 1:ncol(out$scores), sep = "")
  colnames(out$loadings) <- paste("ForeC", 1:ncol(out$loadings), sep = "")
  
  rownames(out$loadings) <- colnames(series)
  class(out$loadings) <- "loadings"
  
  if (ncol(out$loadings) == nrow(out$loadings)) {
    out$mixing <- solve(out$loadings)
  } else {
    out$mixing <- NULL
  }
  
  out$spectrum_method <- spectrum_method
  out$Omega <- omega.scores #Omega(out$scores, spectrum_method = spectrum_method, threshold = threshold)
  out$Omega2[n.comp] = out$Omega[n.comp]

  out$Omega.matrix <- out$loadings %*% diag(out$Omega) %*% t(out$loadings)
  out$threshold <- threshold
  out$sdev <- out$Omega
  out$lambdas <- out$Omega
  class(out) <- c("foreca", "princomp")
  invisible(out)
}


#' @rdname foreca.EM
#' @keywords manip
#' @param series an \eqn{T \times n} array containg a multivariate time series. 
#' Can be a \code{matrix}, \code{data.frame}, or a multivariate \code{ts} object
#' @param f_U multivariate spectrum of class \code{mvspectrum} with 
#' \code{normalize = TRUE}
#' @param spectrum_method string; method to estimate the spectrum. See 
#' \code{\link{mvspectrum}}.
#' @param entropy_method string; method to estimate the entropy. See 
#' \code{\link{discrete_entropy}}.
#' @param tol tolerance level for convergence
#' @param max_iter maximum number of iterations
#' @param kernel an R function \code{kernel = function(lambda) ...} that 
#' weights the frequencies; if \code{NULL} 
#' (default) all frequencies get equal weight.
#' @param nstart how many random starts should be done to avoid local minima?
#' @param threshold set spectral density values below \code{threshold} to 0
#' @param smoothing indicator; if \code{TRUE} then spectrum will be additionally
#' smoothed using \code{\link[mgcv]{gam}}. See \code{\link{mvspectrum}} for details.
#' @param ... other arguments passed to \code{\link{mvspectrum}}
#' @export
#' @examples
#' \dontrun{
#' XX = diff(log(EuStockMarkets)) * 100
#' one_weight = foreca.EM.opt_comp(XX, smoothing = FALSE)
#' 
#' plot(one_weight)
#' }
#' 
foreca.EM.opt_comp <- function(series,
                                       f_U = NULL,
                                       spectrum_method = "multitaper", 
                                       entropy_method = "MLE",
                                       tol = 1e-4, 
                                       max_iter = 100, 
                                       kernel = NULL, 
                                       nstart = 5, 
                                       threshold = 0,
                                       smoothing = FALSE,
                                       ...) {
  if (!is.ts(series)) {
    series <- ts(series)
  }
    
  nseries <- ncol(series)
  PW <- whiten(series)
  UU <- PW$U
  #print(UU[1:3, ])
  ff_UU <- f_U
  if (is.null(ff_UU)){
    ff_UU <- mvspectrum(UU, method = spectrum_method, normalize = TRUE, ...)
    if (!is.null(kernel)) {
      ff_UU <- sweep(ff_UU, 1, kernel(attr(ff_UU, "frequency")), FUN = "*")
      ff_UU <- normalize_mvspectrum(ff_UU)
    }
  }
  
  Omega_best <- 0
  HH_best <- Inf
  WW_best <- NA
  HH_tries <- rep(NA, nstart)
  HH_trace <- NA
  converged <- FALSE
  
  for (TRY in 1:nstart) {
    temp_converged <- FALSE
    
    WW <- matrix(NA, ncol = nseries, nrow = 1)
    if (TRY == 1) {
      WW[1, ] <- PW$Sigma0.5[, 1]
      WW[1, ] <- WW[1, ]/base::norm(WW[1, ], "2")
    } else if (TRY == 2) {
      WW[1, ] <- foreca.EM.init_weightvector(UU, ff_UU, method = "SFA_slow")
    } else if (TRY == 3) {
      WW[1, ] <- foreca.EM.init_weightvector(UU, ff_UU, method = "SFA_fast")
    } else if (TRY == 4){
      WW[1, ] <- foreca.EM.init_weightvector(UU, ff_UU, method = "SFA_slow", 
                                             lag = frequency(series))
    } else if (TRY == 5) {
      WW[1, ] <- foreca.EM.init_weightvector(UU, ff_UU, method = "max")
    } else if (TRY == 6) {
      WW[1, ] <- foreca.EM.init_weightvector(UU, ff_UU, method = "unif")
    } else if (TRY == 7) {
      WW[1, ] <- foreca.EM.init_weightvector(UU, ff_UU, method = "norm")
    } else {
      WW[1, ] <- foreca.EM.init_weightvector(UU, ff_UU, method = "cauchy")
    }
  
    WW[1, ] <- WW[1, ]/sign(WW[1, which.max(abs(WW[1, ]))])
        
    HH <- rep(NA, 1)
    OO <- rep(NA, 1)
    HH_after_M <- rep(NA, 1)
    OO_after_M <- rep(NA, 1)
    
    #print((UU %*% WW[1, ])[1:3])
    
    for (ii in 1:max_iter) {
      f_current <- foreca.EM.E_step(f_U = ff_UU, weights = WW[ii, ])
      #f_current_direct <- mvspectrum(UU %*% WW[ii, ], method = spectrum_method,
      #                               normalize = TRUE, smoothing = smoothing)
      #matplot(cbind(f_current, f_current_direct), log = "y")
      #f_current <- f_current_direct
      #Sys.sleep(0.1)
      HH[ii] <- foreca.EM.h(weights_new = WW[ii, ], f_U = ff_UU, 
                            f_current = f_current, base = NULL)
      
      WW <- rbind(WW, foreca.EM.M_step(ff_UU, f_current, 
                                       minimize = TRUE)$vector)
      
      HH_tries[TRY] <- HH[length(HH)]
      if (ii > 1) {
        rel_improv = HH[ii-1] / HH[ii]
        # cat('Increase in h(w):', round(HH[ii] - HH[ii-1], 4), '\n')
        if (ii > 1 && (abs(rel_improv - 1) < tol)) {
          # convergence criterion
          temp_converged <- TRUE
          break
        }
      }
    }  # end of iter

    Omega_temp <- Omega(UU %*% WW[nrow(WW),], spectrum_method = spectrum_method, 
                        threshold = threshold) #100*(1 - out$h) #Omega_best
    #print(Omega_temp)
    #if (HH[length(HH)] < HH_best) {
    if (Omega_temp > Omega_best) {
      WW_best <- WW
      HH_trace <- HH
      HH_best <- HH[length(HH)]
      Omega_best <- Omega_temp
      converged <- temp_converged
      best_try <- TRY
      first_rel_improv <- HH_best / HH[1]
    }
    if (!temp_converged) {
      TRY <- TRY - 1
    }
  }  # end of TRY
  
  if (!converged) {
    warning("Convergence has not been reached. Please try again.")
  }
  
  out <- list()
  out$n.obs <- nrow(series)
  #out$x <- series
  out$tol <- tol
  out$weights <- WW_best
  rownames(out$weights) = paste("Iter", 1:nrow(out$weights))
  colnames(out$weights) = colnames(series)
  
  out$best.try <- best_try
  
  out$h.trace <- HH_trace
  out$h <- HH_best
  out$h.rel.improvement <- first_rel_improv
  out$best.weights <- as.matrix(out$weights[nrow(out$weights), ])
  
  out$iterations <- nrow(WW_best)
  out$converged <- converged
  
  out$loadings <- PW$Sigma0.5inv %*% out$best.weights
  out$Sigma0.5 <- PW$Sigma0.5
  out$Sigma0.5inv <- PW$Sigma0.5inv
  #out$loadings <- sweep(out$loadings, 2, 
  #                      apply(out$loadings, 2, 
  #                            function(x) base::norm(as.matrix(x), "F")))
  
  out$loadings <- out$loadings * (2 * (out$loadings[1] > 0) - 1)
  
  rownames(out$loadings) <- colnames(series)
  class(out$loadings) <- "loadings"
  
  out$scores <- ts(scale(series, scale = FALSE, center = TRUE) %*% out$loadings)
  
  out$Omega <- Omega_best #Omega(out$scores, spectrum_method = spectrum_method, 
                     #threshold = threshold) #100*(1 - out$h) #Omega_best
  out$h.tries <- HH_tries
  out$spectrum_method <- spectrum_method
  out$entropy_method <- entropy_method
  out$smoothing <- smoothing
  out$best.f <- foreca.EM.E_step(f_U = ff_UU, weights = out$best.weights)
  
  class(out) <- "foreca.EM.opt_comp"
  invisible(out)
} 

