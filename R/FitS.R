
fitS <- function(dataIn, xColIndex=NULL, yColIndex=NULL, slopeIn=NULL, depth=1) {
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

  d <- d[order(d$x, decreasing=FALSE),]
  changePt_low <- min(d$x)
  changePt_high <- max(d$x)
  val_low <- min(d$y)
  val_high <- max(d$y)

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

  #standard error of the difference of means
  retObj$stdErrorDiff <- sqrt(retObj$covMat[1,1] + retObj$covMat[2,2] - 2*retObj$covMat[1,2])

  #Binary Segmentation
  retObj$leftPartition <- NULL
  retObj$leftPartition <- NULL
  
  

  keep_dividing <- TRUE
  # When to stop going
  if (depth <= 1 || nrow(dataIn) <= 3) {
    keep_dividing <- FALSE
  }

  # recursion
  if (keep_dividing) {
    cpIndex <- 3
    if (retObj$slopeGenerated) {
      cpIndex <- 4
    }
    cutPoint <- 0
    while (cutPoint < nrow(d) && d$x[cutPoint+1] < retObj$pars[[cpIndex]]) {
      cutPoint <- cutPoint + 1
    }
    retObj$leftPartition <- fitS(d[1:cutPoint,], 1, 2, slopeIn, depth-1)
    retObj$rightPartition <- fitS(d[(cutPoint+1):nrow(d),], 1, 2, slopeIn,
                                  depth-1)
  }
  
  cp_traverser <- function(current_obj) {
    cps <- c()
    if (!is.null(current_obj$leftPartition)) { #if the current_obj has a leftPartition, traverse through it
      cps <- cp_traverser(current_obj$leftPartition)
    }
    cpIndex <- 3
    if (current_obj$slopeGenerated) {
      cpIndex <- 4
    }
    cps <- append(cps, current_obj$pars[[cpIndex]])
    if (!is.null(current_obj$rightPartition)) {
      cps <- append(cps, cp_traverser(current_obj$rightPartition))
    }
    return(cps)
  }

  
  retObj$cp_list <- cp_traverser(retObj)
  
  
  std_error_traverser <- function(current_obj) {
    std_errors <- c()
    if (!is.null(current_obj$leftPartition)) { #if the current_obj has a leftPartition, traverse through it
      std_errors <- std_error_traverser(current_obj$leftPartition)
    }
    cpIndex <- 3 
    if (current_obj$slopeGenerated) {
      cpIndex <- 4
    }
    std_errors <- append(std_errors, current_obj$stdErrorDiff)
    if (!is.null(current_obj$rightPartition)) {
      std_errors <- append(std_errors, std_error_traverser(current_obj$rightPartition))
    }
    return(std_errors)
  }
  
  retObj$std_error_list <- std_error_traverser(retObj)
  
  

  class(retObj) <- c('fittedS')
  retObj

}

print.fittedS <- function(obj, listAllCp=FALSE)
{
   print('point estimates of the alpha_i')
   print(obj$pars)
   print('covariance matrix')
   print(obj$covMat)
   print('standard error of the difference between pre-changepoint and post-changepoint means')
   print(obj$stdErrorDiff)
   if (listAllCp) { #print all POSSIBLE changepoints
     #moving traverser function for now
     
     print('All changepoints listed as')
     print(cp_list)
     
     #print all the corresponding standard errors
     print('All corresponding standard errors (of pre-mean/post-mean differences) listed as')
     print(std_errors_list)
   }
   

  
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





