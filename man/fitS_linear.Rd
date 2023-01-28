\name{fitS_linear}
\alias{fitS_linear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{S curve fit to gaps}
\description{
  Process the incoming data frame and feed it to fitS
}
\usage{
fitS_linear(dataIn, xColIndex, yColIndex, slopeIn,
          depth, family_wise_error_rate, autoTraverse, fitIntercept)
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
    Forward to fitS
  }
  \item{depth}{
    Forward to fitS
  }
  \item{family_wise_error_rate}{
    Forward to fitS
  }
  \item{autoTraverse}{
    Forward to fitS
  }
  \item{fitIntercept}{
    To decide whether intercepts or slopes will be fit
  }
}
\details{
  This function fits a logistic S-curve to given data frame's slopes.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  The same as fitS.
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
function (x)
{
  }
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
