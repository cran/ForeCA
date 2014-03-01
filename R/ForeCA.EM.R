#' @title EM-like algorithm to estimate optimal ForeCA transformation
#' @name foreca.EM
#' @aliases foreca.EM.opt_weightvector
#' @description 
#' 
#' \code{foreca.EM} estimates the optimal transformation of multivariate time
#' \code{series} to obtain forecastable signals using an EM-like algorithm.
#' 
#' @inheritParams foreca.EM.opt_weightvector
#' @param n.comp number of components to be extracted
#' @param plot logical; if \code{TRUE} a plot of the current optimal 
#' solution \eqn{\mathbf{w}_i^*} will be shown and updated in each iteration
#' of the EM-like algorithm.
#' @return
#' A list with similar output as \code{\link[stats]{princomp}}. Signals are
#' ordered from most to least forecastable. 
#' @export
#' @examples
#' \dontrun{
#' XX <- diff(log(EuStockMarkets[-c(1:50),])) * 100
#' plot(ts(XX))
#' foreca <- foreca.EM(XX[,1:3], n.comp = 3)
#' 
#' summary(foreca)
#' plot(foreca)
#' }

foreca.EM <- function(series, 
                      spectrum.method = 
                        c("multitaper", "direct", "lag window", "wosa", "mvspec"),
                      n.comp = 2, kernel = NULL, plot = TRUE, threshold = 0, 
                      control = list(tol = 1e-6, max.iter = 200, num.starts = 10), 
                      smoothing = FALSE, ...) {
  
  spectrum.method <- match.arg(spectrum.method)
  if (!is.ts(series)) {
    series <- ts(series)
  }
  
  PW <- whiten(series)
  UU <- PW$U
  
  kk <- ncol(UU)
  null.space <- diag(1, kk)
  
  ff.UU <- mvspectrum(UU, method = spectrum.method, normalize = TRUE, ...)
  if (!is.null(kernel)) {
    ff.UU <- sweep(ff.UU, 1, kernel(attr(ff.UU, "frequency")), FUN = "*")
    ff.UU <- normalize_mvspectrum(ff.UU)
  }
  
  out <- list(loadings = matrix(NA, ncol = n.comp, nrow = kk),
              h = rep(NA, n.comp),
              h.tries = matrix(NA, ncol = n.comp, nrow = control$num.starts),
              weights = list(),
              Omega2 = rep(NA, n.comp),
              h.orig = rep(NA, n.comp),
              best.tries = rep(NA, n.comp))
  
  for (ii in seq_len(n.comp)) {
    UU.ortho <- as.matrix(UU) %*% null.space
    ff.UU.ortho <- NULL
    # ff.UU.ortho <- array( t(apply(ff.UU, 1, function(mat, x) t(x) %*% mat
    # %*% x, x = null.space)), dim = c(dim(ff.UU)[1], dim(null.space)[2],
    # dim(null.space)[2] ) ) print(UU.ortho[1:2, ])
    if (ii == kk) {
      weight_ortho <- 1
      out$h[ii] <- spectral_entropy(UU.ortho, spectrum.method = spectrum.method, 
                                    threshold = threshold)
      
    } else {
      foreca.opt.wv <- 
        foreca.EM.opt_weightvector(UU.ortho, f.U = ff.UU.ortho, 
                                   spectrum.method = spectrum.method, 
                                   kernel = kernel, control = control,
                                   threshold = threshold, 
                                   smoothing = smoothing)
      out$best.tries[ii] <- foreca.opt.wv$best.try
      out$h.orig[ii] <- foreca.opt.wv$h.trace[1]
      weight_ortho <- foreca.opt.wv$loadings
      out$weights[[ii]] <- foreca.opt.wv$weights
      out$Omega2[ii] <- foreca.opt.wv$Omega
      
      out$h[ii] <- foreca.opt.wv$h
      out$h.tries[, ii] <- foreca.opt.wv$h.tries
      
      # plot results
      if (plot) {
        plot(foreca.opt.wv, main = paste("Component", ii))
        Sys.sleep(0.25)
      }
    }
    
    out$loadings[, ii] <- null.space %*% weight_ortho
    out$loadings[, ii] <- 
      out$loadings[, ii]/base::norm(as.matrix(out$loadings[, ii]), "F")
    null.space <- Null(out$loadings[, 1:ii])
  }
  
  if (any(is.na(out$loadings))) {
    out$loadings[, kk] <- Null(out$loadings[, -kk])
    out$loadings[, kk] <- 
      out$loadings[, kk]/base::norm(as.matrix(out$loadings[, kk]), "F")
  }
  
  out <- c(out, 
           list(x = series,
                control = control,
                scores2 = ts(as.matrix(UU) %*% out$loadings, 
                             start = start(series), end = end(series)),
                loadings = PW$Sigma0.5inv %*% out$loadings))
  
  out$loadings <- sweep(out$loadings, 2, apply(out$loadings, 2, norm, "2"), "/")
  out$Sigma0.5 <- PW$Sigma0.5
  out$Sigma0.5inv <- PW$Sigma0.5inv
  # make the sign of the loadings consistent for repeated runs
  out$loadings <- sweep(out$loadings, 2, sign(out$loadings[1, ]), "*")
  
  out$n.obs <- nrow(series)
  out$scores <- ts(scale(series, scale = FALSE, center = TRUE) %*% out$loadings, 
                   start = start(series), end = end(series))
  
  omega.scores <- Omega(out$scores, spectrum.method = spectrum.method, 
                        threshold = threshold)
  
  omega_order <- order(omega.scores, decreasing = TRUE)
  
  out$scores <- out$scores[, omega_order]
  omega.scores <- omega.scores[omega_order]
  out$loadings <- out$loadings[, omega_order]
  
  colnames(out$scores) <- paste0("ForeC", seq_len(ncol(out$scores)))
  colnames(out$loadings) <- paste0("ForeC", seq_len(ncol(out$loadings)))
  
  rownames(out$loadings) <- colnames(series)
  class(out$loadings) <- "loadings"
  
  if (ncol(out$loadings) == nrow(out$loadings)) {
    out$mixing <- solve(out$loadings)
  } else {
    out$mixing <- NULL
  }
  
  out <- c(out,
           list(spectrum.method = spectrum.method,
                Omega = omega.scores))
  
  out <- c(out,
           list(Omega.matrix = out$loadings %*% diag(out$Omega) %*% t(out$loadings),
                threshold = threshold,
                sdev = out$Omega,
                lambdas = out$Omega))
  
  class(out) <- c("foreca", "princomp")
  invisible(out)
} 


#' @rdname foreca.EM
#' @keywords manip
#' @param series a \eqn{T \times K} array containing a multivariate time series. 
#' Can be a \code{matrix}, \code{data.frame}, or a multivariate \code{ts} object.
#' @param f.U multivariate spectrum of class \code{mvspectrum} with 
#' \code{normalize = TRUE}
#' @param spectrum.method string; method to estimate the spectrum; see 
#' \code{\link{mvspectrum}}.
#' @param entropy.method string; method to estimate the entropy; see 
#' \code{\link{discrete_entropy}}.
#' @param control a list with control parameters for the iterative EM algorithm:
#' i) \code{tol} sets the tolerance level for convergence; ii) \code{max.iter}
#' sets the maximum number of iterations; and iii) \code{num.starts} determines how 
#' many random starts should be done to avoid local minima (default:\code{5})
#' @param kernel an R function \code{kernel = function(lambda) ...} that 
#' weighs the Fourier frequencies; if \code{NULL} (default) all frequencies get
#' equal weight.
#' @param threshold set estimate of spectral density below \code{threshold} to 0
#' (and renormalize so it adds up to one)
#' @param smoothing logical; if \code{TRUE} the spectral density estimate will be
#' smoothed using \code{\link[mgcv]{gam}}. See \code{\link{mvspectrum}} for
#' further details.
#' @param ... other arguments passed to \code{\link{mvspectrum}}
#' @export
#' @examples
#' \dontrun{
#' XX <- diff(log(EuStockMarkets)) * 100
#' one.weight <- foreca.EM.opt_weightvector(XX, smoothing = FALSE)
#' 
#' plot(one.weight)
#' }
#' 
foreca.EM.opt_weightvector <- 
  function(series, f.U = NULL, 
           spectrum.method = 
             c("multitaper", "direct", "lag window", "wosa", 
               "mvspec"), 
           entropy.method = "MLE", 
           control = list(tol = 1e-6, max.iter = 100, num.starts = 10),
           kernel = NULL, threshold = 0, smoothing = FALSE, ...) {
  
  spectrum.method <- match.arg(spectrum.method)
  if (!is.ts(series)) {
    series <- ts(series)
  }
  
  nseries <- ncol(series)
  PW <- whiten(series)
  UU <- PW$U
  # print(UU[1:3, ])
  ff.UU <- f.U
  if (is.null(ff.UU)) {
    ff.UU <- mvspectrum(UU, method = spectrum.method, normalize = TRUE, 
                        ...)
    if (!is.null(kernel)) {
      ff.UU <- sweep(ff.UU, 1, kernel(attr(ff.UU, "frequency")), FUN = "*")
      ff.UU <- normalize_mvspectrum(ff.UU)
    }
  }
  
  Omega.best <- 0
  HH.best <- Inf
  WW.best <- NA
  HH.tries <- rep(NA, control$num.starts)
  HH.trace <- NA
  converged <- FALSE
  
  for (TRY in seq_len(control$num.starts)) {
    temp.converged <- FALSE
    
    WW <- matrix(NA, ncol = nseries, nrow = 1)
    if (TRY == 1) {
      WW[1, ] <- PW$Sigma0.5[, 1]
      WW[1, ] <- WW[1, ]/base::norm(WW[1, ], "2")
    } else if (TRY == 2) {
      WW[1, ] <- initialize_weightvector(UU, ff.UU, method = "SFA.slow")
    } else if (TRY == 3) {
      WW[1, ] <- initialize_weightvector(UU, ff.UU, method = "SFA.fast")
    } else if (TRY == 4) {
      WW[1, ] <- initialize_weightvector(UU, ff.UU, method = "norm")
    } else if (TRY == 5) {
      WW[1, ] <- initialize_weightvector(UU, ff.UU, method = "max")
    } else if (TRY == 6) {
      WW[1, ] <- initialize_weightvector(UU, ff.UU, method = "unif")
    } else if (TRY == 7) {
      WW[1, ] <- initialize_weightvector(UU, ff.UU, method = "SFA.slow", 
                                             lag = frequency(series))
    } else {
      WW[1, ] <- initialize_weightvector(UU, ff.UU, method = "cauchy")
    }
    
    WW[1, ] <- WW[1, ]/sign(WW[1, which.max(abs(WW[1, ]))])
    
    HH <- rep(NA, 1)  
    for (ii in seq_len(control$max.iter)) {
      if (smoothing) {
        f.current <- mvspectrum(UU %*% WW[ii, ], 
                                method = spectrum.method, 
                                normalize = TRUE, smoothing = TRUE)
      } else {
        f.current <- foreca.EM.E_step(f.U = ff.UU, weightvector = WW[ii, ])
      }
      # f.current_direct <- mvspectrum(UU %*% WW[ii, ], method =
      # spectrum.method, normalize = TRUE, smoothing = smoothing)
      # matplot(cbind(f.current, f.current_direct), log = 'y') f.current <-
      # f.current_direct Sys.sleep(0.1)
      HH[ii] <- foreca.EM.h(weightvector.new = WW[ii, ], f.U = ff.UU, 
                            f.current = f.current, 
                            base = NULL)
      
      WW <- rbind(WW, foreca.EM.M_step(ff.UU, f.current, minimize = TRUE)$vector)
      
      HH.tries[TRY] <- HH[length(HH)]
      if (ii > 1) {
        rel_improv <- HH[ii - 1]/HH[ii]
        # cat('Increase in h(w):', round(HH[ii] - HH[ii-1], 4), '\n')
        if (abs(rel_improv - 1) < control$tol) {
          # convergence criterion
          temp.converged <- TRUE
          break
        }
      }
    }  # end of iter
    
    Omega.tmp <- Omega(UU %*% WW[nrow(WW), ], spectrum.method = spectrum.method, 
                       threshold = threshold)  #100*(1 - out$h) #Omega.best

    # print(Omega.tmp) if (HH[length(HH)] < HH.best) {
    if (Omega.tmp > Omega.best) {
      WW.best <- WW
      HH.trace <- HH
      HH.best <- HH[length(HH)]
      Omega.best <- Omega.tmp
      converged <- temp.converged
      best.try <- TRY
      first.relative.improv <- HH.best/HH[1]
    }
    if (!temp.converged) {
      TRY <- TRY - 1
    }
  }  # end of TRY
  
  if (!converged) {
    warning("Convergence has not been reached. Please try again.")
  }
  
  rownames(WW.best) <- paste("Iter", seq_len(nrow(WW.best)))
  colnames(WW.best) <- colnames(series)
  
  out <- list(n.obs = nrow(series),
              control = control,
              weights = WW.best,
              best.try = best.try,
              h.tries = HH.tries,
              h.trace = HH.trace,
              h = HH.best,
              h.rel.improvement = first.relative.improv,
              iterations = nrow(WW.best),
              converged = converged,
              spectrum.method = spectrum.method,
              entropy.method = entropy.method,
              smoothing = smoothing)
  
  out$best.weights <- as.matrix(out$weights[nrow(out$weights), ])
  out$loadings <- PW$Sigma0.5inv %*% out$best.weights
  out$Sigma0.5 <- PW$Sigma0.5
  out$Sigma0.5inv <- PW$Sigma0.5inv

  # make the sign of the first loading consistent across runs
  out$loadings <- out$loadings * (2 * (out$loadings[1] > 0) - 1)
  rownames(out$loadings) <- colnames(series)
  class(out$loadings) <- "loadings"
  
  out <- c(out,
           list(scores = 
                  ts(scale(series, scale = FALSE, center = TRUE) %*% out$loadings),
                Omega = Omega.best,  #Omega(out$scores, spectrum.method = spectrum.method, threshold = threshold) #100*(1 - out$h) #Omega.best
                best.f = foreca.EM.E_step(f.U = ff.UU, weightvector = out$best.weights)))
  
  class(out) <- "foreca.EM.opt_weightvector"
  invisible(out)
} 
