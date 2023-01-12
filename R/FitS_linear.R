fitS_linear <- function(dataIn, xColIndex=NULL, yColIndex=NULL, slopeIn=NULL,
                     depth=1, family_wise_error_rate = .05, autoTraverse=TRUE
                     ) {
  if (is.ts(dataIn)) {
    d <- data.frame(x=time(dataIn),y=as.numeric(dataIn))
  } else if (is.vector(dataIn)) {
    d <- data.frame(x=1:length(dataIn),y=dataIn)
  } else {
    d <- data.frame(x=dataIn[[xColIndex]], y=dataIn[[yColIndex]])
  }

  d <- d[order(d$x, decreasing=FALSE),]
  cut_d_rear <- d[2:nrow(d),]
  cut_d_front <- d[1:(nrow(d)-1),]
  gap_d_in <- data.frame(x=cut_d_rear$x,
                         b1=(cut_d_rear$y-cut_d_front$y)/(cut_d_rear$x-cut_d_front$x))
  gap_d_in$b0 <- cut_d_rear$y - (gap_d_in$x * gap_d_in$b1)
  return(fitS(dataIn=gap_d_in, xColIndex=1, yColIndex=2, slopeIn=slopeIn,
              depth=depth, family_wise_error_rate=family_wise_error_rate,
              autoTraverse=autoTraverse))
}
