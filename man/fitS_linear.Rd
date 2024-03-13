\name{fitS_linear}
\alias{fitS_linear}
\title{S-Curve Fit, Linear Models }
\description{Estimation of location and magnitudes of change in
intercept and slope, for piecewise linear models.
}
\usage{
fitS_linear(dataIn,xColIndex=NULL,yColIndex=NULL)
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


