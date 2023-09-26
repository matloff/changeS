\name{plot.fittedS}
\alias{plot.fittedS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the graphs}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates a plot (ggplot) of the data passed as argument to the 'fitS' function. Superimposes the fitted S-curve on the data
}
\usage{
%%plot.fittedS(obj)
\method{plot}{fittedS}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
object returned by 'fitS' function

}
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


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


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
