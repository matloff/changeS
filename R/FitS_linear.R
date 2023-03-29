fitS_linear <- function(dataIn, xColIndex=NULL, yColIndex=NULL,
                     family_wise_error_rate = .05) {
  if (is.ts(dataIn)) {
    d <- data.frame(x=time(dataIn),y=as.numeric(dataIn))
  } else if (is.vector(dataIn)) {
    d <- data.frame(x=1:length(dataIn),y=dataIn)
  } else {
    d <- data.frame(x=dataIn[[xColIndex]], y=dataIn[[yColIndex]])
  }

  big_linear_guy <- function(b1, h1, s1, c, b2, h2, s2, x) {
    part1 <- x*(b1+(h1-b1)/(1+exp(-s1*(x-c))))
    part2 <- b2+(h2-b2)/(1+exp(-s2*(x-c)))
    return(part1+part2)
  }

  ret <- nls_multstart(y~big_linear_guy(b1, h1, s1, c, b2, h2, s2, x=x),
                       d, iter=rep(5,7),
                       lower=c(b1=-Inf,h1=-Inf,s1=0,c=min(d$x),b2=-Inf,h2=-Inf,s2=0),
                       start_lower=c(0,0,0,min(d$x),0,0,0),
                       start_upper=c(20,20,10^4,max(d$x),20,20,10^4),
                       supp_errors="Y")

  retObj <- list(nlsOut=ret)  # returned object from nls_multstart
  retObj$pars <- ret$m$getAllPars()  # pre-, post-means, maybe slope, changept
  retObj$d <- d  # x and y
  retObj$d$fitted <- predict(ret) #  generate predictions and add to dataframe 'd'

  class(retObj) <- c("fittedS_linear", "fittedS")
  retObj
}
