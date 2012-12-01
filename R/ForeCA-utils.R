#' @title ForeCA plots, summary, print, etc methods 
#' @name ForeCA-utils
#' @param ... additional arguments passed to 
#' \code{\link[stats]{biplot.princomp}}, \code{\link[stats]{biplot.default}},
#' \code{\link[graphics]{plot}}, or \code{\link[base]{summary}}.
#' @examples
#' # see examples in 'ForeCA.EM'
#'
NULL

#' @rdname ForeCA-utils
#' @method plot ForeCA
#' @description
#' \code{plot.ForeCA} plots the output of a \code{"ForeCA"} object. It is a 
#' visual summary of the analysis with biplots, screeplots, and white noise 
#' tests.
#' @keywords manip
#' @param x an object of class \code{"ForeCA"} or 
#' \code{"ForeCA.EM.one_weightvector"}
#' @export
#' 

plot.ForeCA <- function(x, lag = 10, alpha = 0.05, ...) {
  object <- x
  n.comp <- ncol(object$scores)

  SO <- summary(object, lag = lag, alpha = alpha)
  
  if (all(SO$p.value > alpha)) {
    select_comps <- 1
  } else {
    select_comps <- which(SO$p.value < alpha)
  }
  k_opt <- length(select_comps)
  
  layout(matrix(1:6, ncol = 3, byrow = TRUE))
  par(mar = c(4, 4.5, 3, 2))
  biplot(object)
  
  barplot(as.vector(object$Omega), main = "Forecastability", 
          names.arg = names(object$Omega), 
          ylab = expression(hat(Omega)(x[t]) ~ " (in %)"), 
          ylim = c(0, max(c(SO$Omega.orig, object$Omega)) * 1.05))
  
  abline(h = object$Omega[1], lty = 2, col = 4)
  abline(h = object$Omega[n.comp], lty = 3, col = 4)
  
  barplot(as.vector(SO$Omega.orig), main = "Forecastability", 
          names.arg = names(SO$Omega.orig), 
          ylab = expression(hat(Omega)(x[t]) ~ " (in %)"), 
          ylim = c(0, max(c(SO$Omega.orig, object$Omega)) * 1.05))
  abline(h = object$Omega[1], lty = 2, col = 4)
  abline(h = object$Omega[n.comp], lty = 3, col = 4)
  
  if (ncol(object$scores) >= 4) {
    biplot(object, 3:4)
  } else {
    plot.new()
    #plot(1, type = "n", axes = F, xlab = "", ylab = "")
  }
  barplot(SO$p.value, ylim = c(0, max(alpha * 1.1, SO$p.value)), ylab = "", 
          main = "")
  mtext("p-value \n (H0: white noise)", side = 2, line = 2, cex = 0.8)
  title(paste(sum(SO$p.value > alpha), "white noise"))
  # legend('topleft', c(paste(sum(SO$p.value > alpha), ' white noise')))
  abline(h = alpha, lwd = 2, lty = 2, col = 4)
  
  barplot(SO$p.value.orig, ylim = c(0, max(alpha * 1.1, SO$p.value.orig)), 
          ylab = "", main = "")
  mtext("p-value \n (H0: white noise)", side = 2, line = 2, cex = 0.8)
  title(paste(sum(SO$p.value.orig > alpha), "white noise"))
  abline(h = alpha, lwd = 2, lty = 2, col = 4)
}

#' @rdname ForeCA-utils
#' @method plot ForeCA.EM.one_weightvector
#' @description
#' \code{plot.ForeCA.EM.one_weightvector} gives a visual summary of the EM-like
#' algorithm to obtain a weight-vector \eqn{\mathbf{w}_i^*}.  It gives trace
#' plots of the objective function (\code{\link{ForeCA.EM.h}}) and of the weight
#' vector, and it shows the transformed signal \eqn{y_t^*} 
#' along with its spectral density estimate \eqn{\widehat{f}_y(\omega_j)}.
#' @keywords manip
#' @param main an overall title for the plot: see \code{\link[graphics]{title}}.
#' @export
#' 

plot.ForeCA.EM.one_weightvector <- function(x, main = "", ...) {
  
  object <- x
  max_iter <- object$iterations
  temp_txt <- substitute(list(paste(hat(Omega) == a1, "%")), 
                         list(a1 = round(object$Omega, 2)))
  
  if (object$smoothing){
    spec_temp = c(object$best.f)
    mod_temp = gam(spec_temp~s(c(1:length(spec_temp))))
    object$best.f.smoothed = mod_temp$fit
    object$best.f.smoothed = object$best.f.smoothed / mean(object$best.f.smoothed) / 2
  } else {
    object$best.f.smoothed = NULL
  }
  
  layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = FALSE), heights = c(2, 3))
  par(mar = c(0, 4.5, 2, 1))
  plot(-1 + 1:max_iter, c(object$h.trace, object$h.trace[max_iter - 1]), 
       type = "l", lwd = 2, ylab = "", xlab = "", axes = FALSE, 
       main = main, ...)
  axis(2)
  box()
  grid()
  points(-1 + 1:max_iter, c(object$h.trace, object$h.trace[max_iter - 1]), 
         pch = 19)
  
  # legend('bottomright', paste(expression(Omega)))
  mtext(expression(paste("h(w |", hat(f)[U](omega[j]), ")")), side = 2, 
        line = 2, cex = 1.2)
  #mtext(temp_txt, side = 3, adj = 1, line = -2, cex = 1.5)
  par(mar = c(3, 4.5, 0.5, 1))
  matplot(-1 + 1:max_iter, object$weights, type = "l", ylab = "", xlab = "", 
          axes = FALSE, lwd = 2, ...)
  axis(2)
  axis(1, at = -1 + 1:max_iter)
  mtext("weights", side = 2, line = 2, cex = 1.2)
  mtext("Iteration", side = 1, line = 2, cex = 1.2)
  box()
  grid()
  abline(h = 0)
  matpoints(-1 + 1:max_iter, object$weights, pch = 19, cex = 1.2)
  
  par(mar = c(0, 4.5, 2, 1))
  plot(object$score, ylab = "", xlab = "", axes = FALSE, type = "l")
  box()
  grid()
  axis(2)
  
  par(mar = c(3, 4.5, 0.5, 1))
  plot(seq(0, 0.5, length = length(object$best.f)), object$best.f, type = "l",
       ylab = "", xlab = "", log = "y")
  lines(seq(0, 0.5, length = length(object$best.f)), object$best.f.smoothed,
        col = 4, lwd = 2)
  abline(h = 0.5, col = 2, lty = 2, lwd = 2)
  box()
  grid()
  mtext(expression(hat(f)(omega[j])), side = 2, 
        line = 2, cex = 1.2)
  mtext(temp_txt, side = 3, adj = 1, line = -2, cex = 1.5)

}

#' @rdname ForeCA-utils
#' @method summary ForeCA
#' @description
#' \code{summary.ForeCA} is a \code{\link[base]{summary}} method for a 
#' \code{\link{ForeCA}} output.
#' 
#' @keywords manip
#' @param object an object of class \code{"ForeCA"}
#' @param alpha significance level for testing white noise in \code{\link[stats]{Box.test}}
#' @param lag how many lags to test in \code{\link[stats]{Box.test}}
#' @export
#' 

summary.ForeCA <- function(object, lag = 10, alpha = 0.05, ...) {
  aux_pvalues <- function(series) {
    Box.test(series, lag = lag, type = "Ljung-Box")$p.value
  }
  
  PVALUES <- apply(object$scores, 2, aux_pvalues)
  PVALUES.orig <- apply(object$x, 2, aux_pvalues)
  
  out <- list()
  out$p.value <- round(PVALUES, 4)
  out$p.value.orig <- round(PVALUES.orig, 4)
  out$Omega <- object$Omega
  out$Omega.orig <- Omega(object$x, spectrum_method = object$spectrum_method, 
                          threshold = object$threshold)
  out$selected <- which(out$p.value < alpha)
  out$alpha <- alpha
  out$lag <- lag
  return(out)
}


#' @rdname ForeCA-utils
#' @method biplot ForeCA
#' @description
#' \code{biplot.ForeCA} plots a biplot from the output of \code{\link{ForeCA}} 
#' (a wrapper around \code{\link[stats]{biplot.princomp}}).
#' @keywords hplot
#' @export
#'

biplot.ForeCA <- function(x, ...) {
  object_princomp <- x
  class(object_princomp) <- "princomp"
  biplot(object_princomp, ...)
  abline(h = 0)
  abline(v = 0)
}
