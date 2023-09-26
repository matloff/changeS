\name{fitS}
\alias{fitS}
\title{S-curve fit}
\description{
  Fits an S-curve to given data points.
}
\usage{
fitS(dataIn, xColIndex, yColIndex, slopeIn, depth, family_wise_error_rate, autoTraverse)
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
  \item{slopeIn}{
    A fixed slope for the model
  }
  \item{depth}{
    A bound of depth of binary segmentation
  }
  \item{family_wise_error_rate}{
    Alpha value for determining confidence intervals
  }
  \item{autoTraverse}{
    To decide if binary segmentations automatically go deeper
  }
}
\details{
  This function fits a logistic S-curve to given data frames.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  It's an nls object for now.
}
\references{
  https://cran.r-project.org/web/packages/nls.multstart/nls.multstart.pdf
}
\author{
  Lan Jiang
  Collin Kennedy
  Norm Matloff
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
n <- 500
x <- (1:n)/n
y <- vector(length=n)
trueChangePt <-round(n*2/3)
y[1:trueChangePt] <- rnorm(trueChangePt,10,2)
y[(trueChangePt+1):n] <- rnorm(n-trueChangePt,12.5,2)
d <- data.frame(x=x,y=y)
plot(d)
fitS(d,1,2,10)
fitS(d, 1, 2)


n <- 1000
y <- vector(length = n)
x <- seq(1,n,by = 1)
idx <- c(383,855)
part1 <- runif(n = length(x[1:(idx[1]-1)]), min = 0, max = 4) #mean of 2
part2 <- runif(n = length(x[idx[1]:(idx[2]-1)]), min = 0,max = 10) # mean of 5
part3 <- runif(n = length(x[idx[2]:n]), min = 0, max = 2) #mean of 1

#append these values to y vector
y[1:(idx[1]-1)] <- part1
y[idx[1]:(idx[2]-1)] <- part2
y[idx[2]:n] <- part3

#create a dataframe of the data
df <- data.frame(x = x, y = y)

fitS(df, 1, 2, depth=2, autoTraverse = TRUE)
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
