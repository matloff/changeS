\name{fitS_linear}
\alias{fitS_linear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{S curve fit to gaps}
\description{
  Fits the logistic linear model to the given data points.
}
\usage{
fitS_linear(dataIn, xColIndex, yColIndex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataIn}{
    Incoming data frame
  }
  \item{xColIndex}{
    The column index of interested x in the data frame
  }
  \item{yColIndex}{
    The column index of interested y in the data frame
  }
}
\details{
  This function fits a logistic S-curve to given data frame.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{b1} {
    The pre-changepoint slope.
  }
  \item{h1} {
    The post-changepoint slope.
  }
  \item{s1} {
    The change rate of the slope.
  }
  \item{c} {
    The changepoint.
  }
  \item{b2} {
    An indicator to the pre-changepoint intercept.
  }
  \item{h2} {
    An indicator to the post-changepoint intercept.
  }
  \item{s2} {
    The change rate of the intercept.
  }
}
\references{
  https://cran.r-project.org/web/packages/nls.multstart/nls.multstart.pdf
}
\author{
  Lan Jiang
  Collin Kennedy
  Norm Matloff
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
n <- 200
x <- 1:n + rep(1000, n)
y <- vector(length=n)
trueChangePt <-round(n*2/5)  # 80
y[1:trueChangePt] <- rnorm(trueChangePt,0,10) + seq(from=20, to=60, length.out=trueChangePt)
y[(trueChangePt+1):n] <- rnorm(n-trueChangePt,0,10) + seq(from=20, to=80, length.out=n-trueChangePt)
d <- data.frame(x=x,y=y)
our_temp <- fitS_linear(d, 1, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
