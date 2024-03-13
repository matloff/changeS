\name{fitS_linear}
\alias{fitS_linear}
\title{S-Curve Fit, Linear Models }
\description{Estimation of location and magnitudes of change in
intercept and slope, for piecewise linear models.
}
\usage{
fitS_linear(dataIn,xColIndex=NULL,yColIndex=NULL,
   family_wise_error_rate = 0.05) 
}

\arguments{
  \item{dataIn}{Data frame or equivalent.}
  \item{xColIndex}{Column number of "x" (typically time).}
  \item{yColIndex}{Column number of "y".}
}

\details{
Linear model analog of \code{fitS}.  Note:  May have long run times.
}

\value{

Object of class \code{"fittedS_linear"}, with components:

  \item{b1}{The pre-changepoint slope.}
  \item{h1}{The post-changepoint slope.}
  \item{s1}{The change in slope.}
  \item{c}{The changepoint.}
  \item{b2}{Pre-changepoint intercept.}
  \item{h2}{Post-changepoint intercept.}
  \item{s2}{Change in intercept.}
}

\author{Lan Jiang, Collin Kennedy, Norm Matloff}

\examples{
   \dontrun{
   n <- 200
   x <- 1:n + rep(1000, n)
   y <- vector(length=n)
   trueChangePt <-round(n*2/5)  # 80
   y[1:trueChangePt] <- rnorm(trueChangePt,0,10) + 
      seq(from=20, to=60, length.out=trueChangePt)
   y[(trueChangePt+1):n] <- rnorm(n-trueChangePt,0,10) + 
      seq(from=20, to=80, length.out=n-trueChangePt)
   d <- data.frame(x=x,y=y)
   fitS_linear(d, 1, 2)
   }
}

