
fitS <- function(dataIn, xColIndex=NULL, yColIndex=NULL, slopeIn=NULL, depth=1,
                 family_wise_error_rate = .05, autoTraverse=TRUE) {

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
     dNames <- names(dataIn)[c(xColIndex,yColIndex)]
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
       start_upper=c(val_high, val_high, 25, changePt_high),
       lower=c(postMean=-Inf,preMean=-Inf,slope=0,changePt=-Inf),
       supp_errors="Y")
  }

  retObj <- list(nlsOut=ret)  # returned object from nls_multstart
  retObj$slopeIn <- slopeIn
  retObj$dNames <- if (exists('dNames')) dNames else NULL
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
  # Get real-time CI for diff between preMean and postMean
  real_time_diff <- retObj$pars[[1]] - retObj$pars[[2]] # post-pre
  error_rate <- family_wise_error_rate
  interval_width <- qnorm(1 - (error_rate/2)) * retObj$stdErrorDiff
  real_time_upper <- real_time_diff + interval_width
  real_time_lower <- real_time_diff - interval_width

  if (autoTraverse) {
    if (real_time_lower < 0 && real_time_upper > 0) {
      keep_dividing <- FALSE
    }
  } else {
    print("The real time confidence interval of difference between pre-changepoint mean and post-changepoint mean based on provided error rate is")
    print(c(real_time_lower, real_time_upper))
    print("Keep the results and potentially dig further? Y/N/D (view the used data)")
    while (TRUE) {
      user_input <- NULL
      user_input <- readline()
      if (user_input == "Y" || user_input == "y") {
        break
      } else if (user_input == "N" || user_input == "n") {
        return(NULL)
      } else if (user_input == "D" || user_input == "d") {
        print(retObj$d)
      } else {
        print("Undefined Input.")
      }
      print("Keep dividing? (Y/N/D for current data used)")
    }
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
    retObj$leftPartition <- fitS(d[1:cutPoint,], 1, 2, slopeIn, depth-1, autoTraverse = autoTraverse)
    retObj$rightPartition <- fitS(d[(cutPoint+1):nrow(d),], 1, 2, slopeIn,
                                  depth-1, autoTraverse = autoTraverse)
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


  #need the pre and post means to calculate confidence intervals
  post_mean_traverser <- function(current_obj) {
    post_means <- c()
    if (!is.null(current_obj$leftPartition)) { #if the current_obj has a leftPartition, traverse through it
      post_means <- post_mean_traverser(current_obj$leftPartition)
    }

    post_means <- append(post_means, current_obj$pars[[1]])
    if (!is.null(current_obj$rightPartition)) {
      post_means <- append(post_means, post_mean_traverser(current_obj$rightPartition))
    }
    return(post_means)
  }

  retObj$post_means <- post_mean_traverser(retObj)

  pre_mean_traverser <- function(current_obj) {
    pre_means <- c()
    if (!is.null(current_obj$leftPartition)) { #if the current_obj has a leftPartition, traverse through it
      pre_means <- pre_mean_traverser(current_obj$leftPartition)
    }

    pre_means <- append(pre_means, current_obj$pars[[2]])
    if (!is.null(current_obj$rightPartition)) {
      pre_means <- append(pre_means, pre_mean_traverser(current_obj$rightPartition))
    }
    return(pre_means)
  }

  retObj$pre_means <- pre_mean_traverser(retObj)

  #default to bonferroni correction for now, will look into potentially allowing user to use different methods
  #for controlling family-wise error rate but this should be sufficient for now...
  #also default to .05 family wise error rate for now
  num_comparisons <- length(retObj$pre_means)

  #lower and upper bounds for the confidence intervals
  lower_bounds <- (retObj$post_means - retObj$pre_means) - qnorm(1 - family_wise_error_rate/num_comparisons)*retObj$std_error_list
  upper_bounds <- (retObj$post_means - retObj$pre_means) + qnorm(1 - family_wise_error_rate/num_comparisons)*retObj$std_error_list

  #included them as attributes for retObj but probably will not need to include these in the future
  retObj$lower_bounds <- lower_bounds
  retObj$upper_bounds <- upper_bounds


  retObj$CI_list  <- vector(mode = 'list', length = length(retObj$cp_list))
  for(i in 1:length(retObj$cp_list)){
    retObj$CI_list[[i]] <- list('Possible Changepoint' = retObj$cp_list[i],
                         'Confidence Interval for Corresponding Difference of Post-Changepoint and Pre-Changepoint Means' = str_glue('({lower},{upper})', lower = retObj$lower_bounds[i], upper = retObj$upper_bounds[i]),
                         'Post-Changepoint Mean' = retObj$post_means[i],
                         'Pre-Changepoint Mean' = retObj$pre_means[i],
                         'Std. Error of Difference' = retObj$std_error_list[i])
  }

  # flip the slope if postMean is smaller
  if (is.null(slopeIn) && retObj$pars[[1]] < retObj$pars[[2]]) {
    retObj$pars[[3]] <- -retObj$pars[[3]]
  }


  class(retObj) <- c('fittedS')
  retObj

}

print.fittedS <- function(x,...)
{
   # adapt LJ code with minimal change
   obj <- x
   listAllCp <- FALSE
   print('point estimates of the alpha_i')
   print(obj$pars)
   if (inherits(obj, "fittedS_linear")) {
     return()
   }
   print('covariance matrix')
   print(obj$covMat)
   print('standard error of the difference between pre-changepoint and post-changepoint means')
   print(obj$stdErrorDiff)
   if (listAllCp) { #print all POSSIBLE changepoints
     #moving traverser function for now

     print('All changepoints listed as')
     print(obj$cp_list)

     #print all the corresponding standard errors
     print('All corresponding standard errors (of pre-mean/post-mean differences) listed as')
     print(obj$std_error_list)
   }
}


summary.fittedS <- function(object,...){
  summary(object$nlsOut)
}

plot.fittedS <- function(x,...)
{

   z <- x
   dn <- z$dNames
   xlb <- if (!is.null(dn)) dn[1] else 'x'
   ylb <- if (!is.null(dn)) dn[2] else 'y'
   plot(z$d$x,z$d$y,cex=0.4,xlab=xlb,ylab=ylb)
   prs <- z$pars
   minX <- min(z$d$x)
   maxX <- max(z$d$x)
   changePt <- prs['changePt']
   preMean <- prs['preMean']
   postMean <- prs['postMean']

   if (is.null(z$slopeIn)) {  # gradual change
      graphics::title('Gradual Change')
      a4 <- z$slope
      h <- function(t,preMean,postMean,changePt,a4) 
         preMean + (postMean-preMean) / (1+exp(-a4*(t-changePt))) 
      g <- function(t) h(t,preMean,postMean,changePt,a4) 
      curve(g,minX,maxX,add=T,col='blue')
   } else {  # abrupt change
      graphics::title('Abrupt Change')
      firstFlat <- prs['preMean']
      lines(c(minX,changePt),c(firstFlat,firstFlat))
      secondFlat <- prs['postMean']
      lines(c(changePt,maxX),c(secondFlat,secondFlat))
   }

}


###  plot.fittedS <- function(x,...){
 
###    # adapt LJ code with minimal change
###    obj <- x
###    title <- "Changepoint Plot"
###    #should return a ggplot object
###    cpIndex <- 3
###    if (obj$slopeGenerated) {
###      cpIndex <- 4
###    }
###  
###    annot <- data.frame(
###      x = c(obj$par[[cpIndex]],(obj$par[[cpIndex]]+min(obj$d$x))/2, 
###         (obj$par[[cpIndex]]+max(obj$d$x))/2), #cp,pre, post
###      y = c(max(obj$d$y), obj$pars[[2]]/2, (max(obj$d$y) - obj$pars[[1]])/2 ),
###      label = c(str_glue("Estimated CP: {cp}", 
###         cp = round(obj$par[[cpIndex]],2)),  #changepoint
###                str_glue("Pre-CP Mean = {pre}", 
###                pre = round(obj$pars[[2]],3)),     #pre-cp mean
###                str_glue("Post-CP Mean = {post}", 
###                post = round(obj$pars[[1]],3)))  #post-cp mean
###    )
###  
###    x <- annot$x; y <- annot$y; label <- annot$label
###    ggplot(data = obj$d, mapping = aes(x = x, y = y))+
###      geom_point(alpha = .8, color = 'black')+
###      geom_line(data = obj$d, mapping = aes(x = x, y = fitted, color = "S-fit"))+
###      ggtitle(title)+
###      scale_color_manual(name = "Curve Type", values = "blue")+
###      geom_vline(xintercept = obj$par[[cpIndex]], color = "lightblue", linetype = "dashed")+
###      theme_minimal()+
###      geom_label(data = annot, aes(x = x, y = y, label = label), color = "orange", fontface = "bold")+
###      theme(plot.title = element_text(face = "bold", size = 15, hjust= .5),
###            plot.subtitle = element_text(hjust = .5))
###  
###  }





