#' @title Implementation of Forecastable Component Analysis (ForeCA)
#'
#' @description
#' Forecastable Component Analysis (ForeCA) is a novel dimension reduction (DR) 
#' technique for multivariate time series.  ForeCA finds a linar combination
#'  \eqn{y_t = \mathbf{w}' \mathbf{X}_t} that is easy to forecast. The measure of 
#' forecastability \eqn{\Omega(x_t)} (\code{\link{Omega}}) is based on the entropy 
#' of the spectral density \eqn{f_y(\lambda)} of \eqn{y_t}: higher entropy means 
#' less forecastable, lower entropy is more forecastable.
#'
#' The main function \code{\link{foreca}} runs ForeCA on a multivariate time 
#' series to find the most forecastable signals.
#'  
#' Even though the current version has most functionality of this R package, some
#' function name conventions might change in future versions. Please consult the
#' \code{NEWS} file for a list of changes.
#'
#'@name ForeCA-package
#'@aliases ForeCA-package ForeCA
#'@docType package
#'@author Author and maintainer: Georg M. Goerg <im@@gmge.org>
#'@references Goerg, G. M. (2013). \dQuote{Forecastable Component
#' Analysis}. Journal of Machine Learning Research (JMLR) W&CP 28 (2): 64-72, 2013.
#' Available at \url{jmlr.org/proceedings/papers/v28/goerg13.html}.
#'@keywords package
#'@examples
#'XX <- ts(diff(log(EuStockMarkets))[-c(1:1000),])
#'Omega(XX)
#'
#'plot(log(lynx,10))
#'Omega(log(lynx,10), spectrum.method = "wosa")
#'
#'\dontrun{
#'ff <- foreca(XX, n.comp = 4, spectrum.method = "wosa")
#'plot(ff)
#'summary(ff)
#'}
#'
NULL



