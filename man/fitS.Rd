\name{fitS}
\alias{fitS}
\title{S-curve fit}
\description{
  Fits an S-curve to given data points.
}
\usage{
fitS(dataIn, xColIndex, yColIndex, slopeIn)
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
library(changepointS)
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
