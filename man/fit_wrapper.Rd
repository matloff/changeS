\name{fit_wrapper}
\alias{fit_wrapper}
\title{S-curve fit}
\description{
  Fits an S-curve to given data points.
}
\usage{
fit_wrapper(data_in, x_col_index, y_col_index, slope_in=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_in}{
    Incoming data frame
  }
  \item{x_col_index}{
    The column index of interested x in the data frame
  }
  \item{y_col_index}{
    The column index of interested y in the data frame
  }
  \item{slope_in}{
    A fixed slope for the model
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
  Collin Kennedy
  Lan Jiang
  Norm Matloff
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
library(changepointS)
n = 10000;
x = (1:n)/n;
y = vector(length=n);
trueChangePt = round(n*2/3);
y[1:trueChangePt] = runif(trueChangePt, 0, 0.6);
y[(trueChangePt+1):n] = runif(n-trueChangePt, 0.3, 1);

d = data.frame(y=y, x=x);

fit_wrapper(d, 2, 1)

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
