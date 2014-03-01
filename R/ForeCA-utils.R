#' @title Plot, summary, and print methods for class 'foreca'
#' @name foreca-utils
#' @aliases plot.foreca summary.foreca biplot.foreca
#' @param ... additional arguments passed to 
#' \code{\link[stats]{biplot.princomp}}, \code{\link[stats]{biplot.default}},
#' \code{\link[graphics]{plot}}, or \code{\link[base]{summary}}.
#' @examples
#' # see examples in 'foreca.EM'
#'
NULL

#' @rdname foreca-utils
#' @method plot foreca
#' @description
#' \code{plot.foreca} shows a visual summary of the ForeCA results with biplots, 
#' screeplots, and white noise tests.
#' @keywords manip
#' @param x an object of class \code{"foreca"} or \code{"foreca.EM.opt_weightvector"}
#' @export
#' 

plot.foreca <- function(x, lag = 10, alpha = 0.05, ...) {
  object <- x
  n.comp <- ncol(object$scores)

  SO <- summary(object, lag = lag, alpha = alpha)
    
  par(mar = c(4, 4.5, 2.5, 2.5))
  layout(matrix(1:6, ncol = 3, byrow = TRUE))
  biplot(object)
  
  # replace 'Series' with 'ForeC'
  names(object$Omega) <- gsub("Series ", "ForeC", names(object$Omega))
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

#' @rdname foreca-utils
#' @method plot foreca.EM.opt_weightvector
#' @description
#' \code{plot.foreca.EM.opt_weightvector} shows how the EM-like algorithm
#' arrived at the i-th optimal a weight-vector \eqn{\mathbf{w}_i^*}.  It 
#' shows trace plots of the objective function (\code{\link{foreca.EM.h}}) 
#' and of the weight vector, and the transformed signal \eqn{y_t^*} along with 
#' its spectral density estimate \eqn{\widehat{f}_y(\omega_j)}.
#' @keywords manip
#' @param main an overall title for the plot: see \code{\link[graphics]{title}}.
#' @param cex.lab size of the axes labels
#' @export
#' 

plot.foreca.EM.opt_weightvector <- function(x, main = "", cex.lab = 1.1, ...) {
  
  object <- x
  total.iter <- object$iterations
  temp.txt <- substitute(list(paste(hat(Omega) == a1, "%")), 
                         list(a1 = round(object$Omega, 2)))
  
  if (object$smoothing) {
    spec.tmp <- c(object$best.f)
    mod.tmp <- gam(spec.tmp ~ s(seq_along(spec.tmp)))
    object$best.f.smoothed <- mod.tmp$fit
    object$best.f.smoothed <- object$best.f.smoothed/mean(object$best.f.smoothed)/2
  } else {
    object$best.f.smoothed <- NULL
  }
  
  par(mar = c(0, 4.5, 2, 0.5))
  layout(matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE), heights = c(2, 3))

  plot(seq(0, total.iter - 1, by = 1), 
       c(object$h.trace, object$h.trace[total.iter - 1]), type = "l", lwd = 2, 
       ylab = "", xlab = "", axes = FALSE, main = main, ...)
  axis(2)
  box()
  grid()
  points(seq(0, total.iter - 1, by = 1), 
         c(object$h.trace, object$h.trace[total.iter - 1]), pch = 19)
  
  mtext(expression(paste("h(w|", hat(f)[U](omega[j]), ")")), side = 2, line = 2, 
        cex = cex.lab)

  plot(object$score, ylab = "", xlab = "", axes = FALSE, type = "l")
  box()
  grid()
  axis(2)
  
  par(mar = c(3.5, 4.5, 1, 0.5))
  matplot(-1 + seq_len(total.iter), object$weights, type = "l", 
          ylab = "", xlab = "", axes = FALSE, lwd = 2, ...)
  axis(2)
  axis(1, at = -1 + seq_len(total.iter))
  mtext("weights", side = 2, line = 2, cex = cex.lab)
  mtext("Iteration", side = 1, line = 2.5, cex = cex.lab)
  box()
  grid()
  abline(h = 0)
  matpoints(-1 + seq_len(total.iter), object$weights, pch = 19, cex = cex.lab)
  
  plot(seq(0, 0.5, length = length(object$best.f)), object$best.f, type = "l", 
       ylab = "", xlab = "", log = "y")
  lines(seq(0, 0.5, length = length(object$best.f)), object$best.f.smoothed, 
        col = 4, lwd = 2)
  abline(h = 0.5, col = 2, lty = 2, lwd = 2)
  box()
  grid()
  mtext(expression(paste("Frequency / 2", pi)), side = 1, line = 2.5, cex = cex.lab)
  mtext(expression(paste(hat(f)(omega[j]), " (log scale)")), side = 2, 
        line = 2, cex = cex.lab)
  mtext(temp.txt, side = 3, adj = 1, line = -2, cex = cex.lab * 1.1)
} 

#' @rdname foreca-utils
#' @method summary foreca
#' @description
#' \code{summary.foreca} computes summary statistics of the ForeCA results.
#' 
#' @keywords manip
#' @param object an object of class \code{"foreca"}
#' @param alpha significance level for testing white noise in 
#' \code{\link[stats]{Box.test}}; default: \eqn{5\%}.
#' @param lag how many lags to test in \code{\link[stats]{Box.test}}; default
#' \eqn{10} lags.
#' @export
#' 
summary.foreca <- function(object, lag = 10, alpha = 0.05, ...) {
  aux_pvalues <- function(series) {
    Box.test(series, lag = lag, type = "Ljung-Box")$p.value
  }
  
  pvals <- apply(object$scores, 2, aux_pvalues)
  pvals.orig <- apply(object$x, 2, aux_pvalues)
  
  out <- list(p.value = round(pvals, 4), 
              p.value.orig = round(pvals.orig, 4), 
              Omega = object$Omega, 
              Omega.orig = 
                Omega(object$x, spectrum.method = object$spectrum.method,
                      threshold = object$threshold), 
              alpha = alpha,
              lag = lag)
  
  out$selected <- which(out$p.value < alpha)
  
  return(out)
} 

#' @rdname foreca-utils
#' @method biplot foreca
#' @description
#' \code{biplot.foreca} shows a biplot of the ForeCA weightvectors
#' (a wrapper around \code{\link[stats]{biplot.princomp}}).
#' @keywords hplot
#' @export
#'

biplot.foreca <- function(x, ...) {
  object.princomp <- x
  class(object.princomp) <- "princomp"
  biplot(object.princomp, ...)
  abline(h = 0)
  abline(v = 0)
}
