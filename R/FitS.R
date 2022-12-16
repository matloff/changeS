
fitS <- function(dataIn, xColIndex=NULL, yColIndex=NULL, slopeIn=NULL) {

  # Preliminary work (should be a smarter way to do this but haven't found yet)
  # library(nls.multstart)
  fitGenLogit_unfixed <- function(topLevel, bottomLevel, slope, changePt, x) {
    return((topLevel-bottomLevel)/(1+exp(-slope*(x-changePt))) + bottomLevel)
  }
  fitGenLogit_fixed <- function(topLevel, bottomLevel, changePt, x) {
    return((topLevel-bottomLevel)/(1+exp(-slopeIn*(x-changePt))) + bottomLevel)
  }

  # dataIn could be: a data frame; a numeric vector; a time series
  # vector
  
  if (is.ts(dataIn)) {
     d <- data.frame(x=time(dataIn),y=as.numeric(dataIn))
  } else if (is.vector(dataIn)) {
     d <- data.frame(x=1:length(dataIn),y=dataIn)
  } else {
     d <- data.frame(x=dataIn[[xColIndex]], y=dataIn[[yColIndex]])
  }

  changePt_low = min(d$x)
  changePt_high = max(d$x)
  val_low = min(d$y)
  val_high = max(d$y)

  # fit
  ret <- NULL
  if (!is.null(slopeIn)) {
    ret <- nls_multstart(y~fitGenLogit_fixed(topLevel, bottomLevel, changePt, 
       x=x), d, iter=c(5,5,5), start_lower=c(val_low, val_low, changePt_low),
       start_upper=c(val_high, val_high, changePt_high), supp_errors="Y")
  } else {
    ret <- nls_multstart(y~fitGenLogit_unfixed(topLevel, bottomLevel, slope, 
       changePt, x=x),
       d, iter=c(5,5,5,5), start_lower=c(val_low, val_low, 0, changePt_low),
       start_upper=c(val_high, val_high, 25, changePt_high), supp_errors="Y")
  }

  retObj <- list(nlsOut=ret)
  retObj$pars <- ret$m$getAllPars()
  retObj$covMat <- vcov(ret)

  sm <- summary(ret)
  retObj$covMat <- sm$sigma^2 * sm$cov.unscaled 

  class(retObj) <- c('fittedS')
  retObj

}

print.fittedS <- function(obj) 
{
   print(obj$pars)
}

