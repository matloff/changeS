
fitS <- function(dataIn, xColIndex=NULL, yColIndex=NULL, slopeIn=NULL) {

  # Preliminary work (should be a smarter way to do this but haven't found yet)
  # library(nls.multstart)
  fitGenLogit_unfixed <- function(postMean, preMean, slope, changePt, x) {
    return((postMean-preMean)/(1+exp(-slope*(x-changePt))) + preMean)
  }
  fitGenLogit_fixed <- function(postMean, preMean, changePt, x) {
    return((postMean-preMean)/(1+exp(-slopeIn*(x-changePt))) + preMean)
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
    ret <- nls_multstart(y~fitGenLogit_fixed(postMean, preMean, changePt,
       x=x), d, iter=c(5,5,5), start_lower=c(val_low, val_low, changePt_low),
       start_upper=c(val_high, val_high, changePt_high), supp_errors="Y")
  } else {
    ret <- nls_multstart(y~fitGenLogit_unfixed(postMean, preMean, slope,
       changePt, x=x),
       d, iter=c(5,5,5,5), start_lower=c(val_low, val_low, 0, changePt_low),
       start_upper=c(val_high, val_high, 25, changePt_high), supp_errors="Y")
  }

  retObj <- list(nlsOut=ret)  # returned object from nls_multstart
  retObj$pars <- ret$m$getAllPars()  # pre-, post-means, maybe slope, changept
  retObj$d <- d  # x and y
  retObj$d$fitted <- predict(ret) #  generate predictions and add to dataframe 'd'

  # covariance matrix of the par estimates
  sm <- summary(ret)
  retObj$covMat <- sm$sigma^2 * sm$cov.unscaled
  retObj$slopeGenerated <- is.null(slopeIn)

  class(retObj) <- c('fittedS')
  retObj

}

print.fittedS <- function(obj)
{
   print('point estimates of the alpha_i')
   print(obj$pars)
   print('covariance matrix')
   print(obj$covMat)
}

summary.fittedS <- function(obj){
  summary(obj$nlsOut)
}

plot.fittedS <- function(obj){
  #should return a ggplot object
  cpIndex <- 3
  if (obj$slopeGenerated) {
    cpIndex <- 4
  }
  ggplot(data = obj$d, mapping = aes(x = x, y = y))+
    geom_point(alpha = .9, color = 'black')+
    geom_line(data = obj$d, mapping = aes(x = x, y = fitted, color = "S-fit"))+
    ggtitle("Changepoint Plot")+
    scale_color_manual(name = "Curve Type", values = "blue")+
    geom_vline(xintercept = obj$par[[cpIndex]], color = "lightblue", linetype = "dashed")+
    theme_minimal()+
    annotate("text", x = obj$par[[cpIndex]], y = max(obj$d$y), label = str_glue("Estimated CP: {cp}", cp = round(obj$par[[cpIndex]],2)))+
    # annotate("text", x = obj$par[[cpIndex]]+20, y = max(y)-150, label = str_glue('(std. error: {std.error})', std.error = round(std_error_cp,3)))+
    # scale_x_continuous(breaks = seq(min(obj$d$x),max(obj$d$x), by = 10))+
    theme(plot.title = element_text(face = "bold", size = 15, hjust= .5),
          plot.subtitle = element_text(hjust = .5))

}

