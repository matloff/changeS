
fitS_linear <- function(dataIn, xColIndex=NULL, yColIndex=NULL,
    plotTitle='') {

  if (is.ts(dataIn)) {
    d <- data.frame(x=time(dataIn),y=as.numeric(dataIn))
  } else if (is.vector(dataIn)) {
    d <- data.frame(x=1:length(dataIn),y=dataIn)
  } else {
    dNames <- names(dataIn)[c(xColIndex,yColIndex)]
    d <- data.frame(x=dataIn[[xColIndex]], y=dataIn[[yColIndex]])
  }

  big_linear_guy <- function(b1, h1, c, b2, h2, x) {
    part1 <- x*(b1+(h1-b1)/(1+exp(-10*(x-c))))
    part2 <- b2+(h2-b2)/(1+exp(-10*(x-c)))
    return(part1+part2)
  }

  ret <- nls_multstart(y~big_linear_guy(b1, h1, c, b2, h2, x=x),
                       d, iter=rep(5,5),
                       lower=c(b1=-Inf,h1=-Inf,c=min(d$x),b2=-Inf,h2=-Inf),
                       start_lower=c(0,0,min(d$x),0,0),
                       start_upper=c(20,20,max(d$x),20,20),
                       supp_errors="Y")

  retObj <- list(nlsOut=ret)  # returned object from nls_multstart
  retObj$pars <- ret$m$getAllPars()  # pre-, post-means, maybe slope, changept
  retObj$d <- d  # x and y
  retObj$dNames <- if (exists('dNames')) dNames else NULL
  retObj$plotTitle <- plotTitle
  retObj$d$fitted <- predict(ret) #  generate predictions and add to dataframe 'd'

  class(retObj) <- c("fittedS_linear", "fittedS")
  retObj
}

plot.fittedS_linear <- function(x,...)
{
   z <- x
   dn <- z$dNames
   xlb <- if (!is.null(dn)) dn[1] else 'x'
   ylb <- if (!is.null(dn)) dn[2] else 'y'
   plot(z$d$x,z$d$y,cex=0.4,xlab=xlb,ylab=ylb)
   # graphics::title('Piecewise Linear Model')
   prs <- z$pars
   minX <- min(z$d$x)
   maxX <- max(z$d$x)
   endPtsX <- c(minX,prs['c'])
   endPtsY <- c(minX,prs['c'])
   lineFtn <- function(t) prs['b2'] + prs['b1']*t
   endPtsY <- lineFtn(endPtsX)
   graphics::lines(endPtsX,endPtsY,col='blue')
   endPtsX <- c(prs['c'],maxX)
   lineFtn <- function(t) prs['h2'] + prs['h1']*t
   endPtsY <- lineFtn(endPtsX)
   graphics::lines(endPtsX,endPtsY,col='red')
   if (z$plotTitle != '') graphics::title(z$plotTitle)

}

