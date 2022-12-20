\name{plot.fittedS}
\alias{plot.fittedS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates a plot (ggplot) of the data passed as argument to the 'fitS' function. Superimposes the fitted S-curve on the data
}
\usage{
plot.fittedS(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
object returned by 'fitS' function

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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Lan Jiang
Collin Kennedy
Norm Matloff
}
\note{
%%  ~~further notes~~
You can modify the output of this function like typical a ggplot object (i.e., modify title, change x and y axis labels, etc)
https://ggplot2.tidyverse.org/reference/ggplot.html
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

library(changepointS)
library(datasets)
library(broom)
library(ggplot2)
library(stringr)
data("Nile")
nile_df <- tidy(Nile) #convert to dataframe from time series

#test function
out <- fitS(nile_df, xColIndex = 1, yColIndex = 2)
plot(out)


## The function is currently defined as
function (obj) 
{
    ggplot(data = out$d, mapping = aes(x = x, y = y)) + geom_point(alpha = 0.9, 
        color = "black") + geom_line(data = out$d, mapping = aes(x = x, 
        y = fitted, color = "S-fit")) + ggtitle("Changepoint Plot") + 
        scale_color_manual(name = "Curve Type", values = "blue") + 
        geom_vline(xintercept = out$par[[4]], color = "lightblue", 
            linetype = "dashed") + theme_minimal() + annotate("text", 
        x = out$par[[4]], y = max(out$d$y), label = str_glue("Estimated CP: {cp}", 
            cp = round(out$par[[4]], 2))) + theme(plot.title = element_text(face = "bold", 
        size = 15, hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
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
