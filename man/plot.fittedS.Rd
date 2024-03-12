\name{plot.fittedS}
\alias{plot.fittedS}

\title{Graphical Analysis of S-Curve fits}

\description{

Generates a plot of the data passed to \code{fiitS}, with 
the fitted S-curve superimposed 

}
\usage{
   plot(fitSoutputobject)
}

\author{ Lan Jiang, Collin Kennedy, Norm Matloff}

\examples{
   \dontrun{
   nileOut <- fitS(nile,1,2,10)
   plot(nileOut)
   }
}

