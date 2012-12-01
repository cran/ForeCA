#' @title Performs Forecastable Component Analysis
#' @name foreca
#' @description 
#' \code{foreca} performs Forecastable Component Analysis on a given 
#' multivariate and returns an object of class \code{foreca}.
#' 
#' @param series multivariate time series
#' @param method what method should be used; currently only the EM-like 
#' algorithm from the original paper is available (see References): 
#' \code{method = "EM"} calls \code{\link{foreca.EM}}.
#' @param ... additional arguments passed to the individual methods (currently
#' only \code{\link{foreca.EM}}).
#' @return
#' A list with similar output as \code{\link[stats]{princomp}}. Signals are
#' ordered from most to least forecastable. 
#' @export
#' @examples
#' \dontrun{
#' XX = diff(log(EuStockMarkets[-c(1:50),])) * 100
#' plot(ts(XX))
#' ff = foreca(XX[,1:3], n.comp = 3)
#' 
#' summary(ff)
#' plot(ff)
#' }

foreca <- function(series, method = "EM", ...) {

  if (method == "EM"){
    out <- foreca.EM(series, ...)
  } else {
    stop(paste("There is no such method as '", method, "' implemented."))
  }
  
  class(out) <- c("foreca", class(out))
  
  invisible(out)
}