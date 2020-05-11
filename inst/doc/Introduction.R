## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(cache = TRUE)

## ----load-packages, cache=FALSE, message = FALSE------------------------------
library(ForeCA)

## ----cite-ForeCA--------------------------------------------------------------
citation("ForeCA")

## ----set-controls-------------------------------------------------------------
# spectrum control
sc <- list(method = "wosa")
# entropy control
ec <- list(prior.weight = 1e-2)

## ----load-eu-stock-markets----------------------------------------------------
data("EuStockMarkets")

## ----eu-stock-markets---------------------------------------------------------
# log-returns in %
ret <- ts(diff(log(EuStockMarkets)) * 100)

## ----plot-returns, cache = FALSE, echo = FALSE--------------------------------
plot(ret)

## ----cor-matrix---------------------------------------------------------------
cor.ret <- cor(ret)
kable(format(cor.ret, digits = 2), 
      caption = "Correlation matrix")
kable(format(solve(cor.ret), digits = 2),
      caption = "Conditional covariance given other variables")

## ----acf-spectra--------------------------------------------------------------
ret.spec <- mvspectrum(ret, spectrum.control = sc)
plot(ret.spec)
layout(matrix(seq_len(ncol(ret)), ncol = 2))
for (nn in colnames(ret)) {
  pacf(ret[, nn], main = nn)
}

## ----omega-eu-stocks----------------------------------------------------------
ret.omega <- Omega(ret, spectrum.control = sc, entropy.control = ec)
ret.omega

## ----foreca-returns-----------------------------------------------------------
mod.foreca.ret <- foreca(ret, n.comp = ncol(ret), 
                         spectrum.control = sc,
                         entropy.control = ec)

## ----show-results-------------------------------------------------------------
mod.foreca.ret  # this uses the print.foreca method

## ----str-foreca-ret, eval = FALSE---------------------------------------------
#  str(mod.foreca.ret) # too much to display; try it out!

## ----foreca-ret-omega---------------------------------------------------------
mod.foreca.ret$Omega

## ----foreca-ret-loadings------------------------------------------------------
mod.foreca.ret$loadings

## ----foreca-returns-summary---------------------------------------------------
plot(mod.foreca.ret)
plot(ts(mod.foreca.ret$scores))

## ----foreca-unit-variance-check-----------------------------------------------
round(colMeans(mod.foreca.ret$scores), 3)
kable(round(cov(mod.foreca.ret$scores), digits = 3))

## ----session-info-------------------------------------------------------------
sessionInfo()

