\name{fitS}
\alias{fitS}

\title{S-Curve Fit}

\description{Finds change(s) in mean.}
\usage{fitS(dataIn,xColIndex=NULL,yColIndex=NULL,slopeIn=NULL,depth=1,
   family_wise_error_rate=0.05,autoTraverse=TRUE)
}

\arguments{

  \item{dataIn}{Data frame or equivalent.}
  
  \item{xColIndex}{Column number of "x" (typically time).}
  
  \item{yColIndex}{Column number of "y".}
  
  \item{slopeIn}{A fixed slope value for the model. Should be set to a
  large value for the abrupt-change case.  If NULL, the algorithm will
  estimate the slope (gradual change case) }
  
  \item{depth}{Upper bound for the depth of the binary segmentation ;
  if this is 1, the algorithm models the situation of (at most) 1
  changepoint.}
  
  \item{family_wise_error_rate}{Nominal alpha value for determining 
  whether to proceed with the binary segmentations to the next split.
  Meaningful only if \code{depth} larger than 1.}
  
  \item{autoTraverse}{If TRUE, do automatic binary segmentation.
  Meaningful only of \code{depth} larger than 1.}

}

\details{Changepoint detection/estimation for changes in mean,
performed by using an S-curve (logistic function) to 
approximate a step function.  This enables asymptotic standard errors,
and associated confidence intervals and tests for changepoint locations
and change magnitudes.  (However, in the multi-changepoint case, the
alpha levels are only nominal.)

Note:  The location of a changepoint is considered to be a continuous
numeric quantity, in contrast to packages such as \code{changepoints}
where the location is integer-valued.
}

\value{A \code{'fittedS'} object, containing estimates, standard errors
and so on.}   

\author{Lan Jiang, Collin Kennedy, Norm Matloff}

\examples{

# real data

# type ?Nile for background information
nile <- data.frame(t=1871:1970, ht=Nile)
fitS(nile,1,2,10)  # abrupt change model

# type ?cancerRates for background information
data(cancerRates)
fitS(cancerRates,1,2)  # gradual change model

# simulated data, changepoint at i = 367
n <- 500
x <- (1:n)/n
y <- vector(length=n)
trueChangePt <-round(n*2/3)
y[1:trueChangePt] <- rnorm(trueChangePt,10,2)
y[(trueChangePt+1):n] <- rnorm(n-trueChangePt,12.5,2)
d <- data.frame(x=x,y=y)
plot(d)
fitS(d,1,2,10) # abrupt
fitS(d, 1, 2) # gradual

\dontrun{

# simulated data, changepoints at  i= 383, 855
n <- 1000
y <- vector(length = n)
x <- seq(1,n,by = 1)
idx <- c(383,855)
part1 <- runif(n = length(x[1:(idx[1]-1)]), min = 0, max = 4) #mean of 2
part2 <- runif(n = length(x[idx[1]:(idx[2]-1)]), min = 0,max = 10) # mean of 5
part3 <- runif(n = length(x[idx[2]:n]), min = 0, max = 2) #mean of 1
y[1:(idx[1]-1)] <- part1
y[idx[1]:(idx[2]-1)] <- part2
y[idx[2]:n] <- part3
df <- data.frame(x = x, y = y)
fitS(df, 1, 2, depth=2, autoTraverse = TRUE)

}

}

