fit_wrapper <- function(data_in, x_col_index, y_col_index, slope_in=NULL) {

  # Preliminary work (should be a smarter way to do this but haven't found yet)
  # library(nls.multstart)
  fitGenLogit_unfixed <- function(topLevel, bottomLevel, slope, changePt, x) {
    return((topLevel-bottomLevel)/(1+exp(-slope*(x-changePt))) + bottomLevel)
  }
  fitGenLogit_fixed <- function(topLevel, bottomLevel, changePt, x) {
    return((topLevel-bottomLevel)/(1+exp(-slope_in*(x-changePt))) + bottomLevel)
  }

  # build the used data_frame
  d <- data.frame(x=data_in[[x_col_index]], y=data_in[[y_col_index]])
  changePt_low = min(d$x)
  changePt_high = max(d$x)
  val_low = min(d$y)
  val_high = max(d$y)

  # fit
  ret <- NULL
  if (!is.null(slope_in)) {
    ret <- nls_multstart(y~fitGenLogit_fixed(topLevel, bottomLevel, changePt, x=x),
                         d, iter=c(5,5,5), start_lower=c(val_low, val_low, changePt_low),
                         start_upper=c(val_high, val_high, changePt_high), supp_errors="Y")
  } else {
    ret <- nls_multstart(y~fitGenLogit_unfixed(topLevel, bottomLevel, slope, changePt, x=x),
                         d, iter=c(5,5,5,5), start_lower=c(val_low, val_low, 0, changePt_low),
                         start_upper=c(val_high, val_high, 25, changePt_high), supp_errors="Y")
  }
  return(ret) # temporarily returning an nls obj
}


