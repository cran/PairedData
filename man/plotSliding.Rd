\name{plotSliding}
\alias{plotSliding}
\title{Sliding square plot}
\description{
Produce a sliding square plot for paired data which mixes the usual scatterplot 
with the tukey mean-difference plot.
}
\usage{
plotSliding(df,condition1,condition2,xlab = "", ylab = "", ...)
}

\arguments{
  \item{df}{Data frame}
  \item{condition1}{name of the variable corresponding to measurements in the "first" condition}
 \item{condition2}{name of the variable corresponding to measurements in the "second" condition}
 \item{xlab}{new label for X axis}
\item{ylab}{new label for Y axis}
 \item{...}{Arguments to be passed to methods, such as graphical parameters (see par)}
}


\references{
\itemize{
\item P.R. Rosenbaum. Exploratory plot for paired data. American Statistician, 43, 108-110, 1989.
\item J.S. Pontius and R.M. Schantz. Graphical analyses of a twoperiod crossover design. The American Statistician, 48, 249-253, 1994.
\item R.M. Pruzek and J.E. Helmreich. Enhancing dependent sample analyses with graphics. Journal of Statistics Education, 17, 2009.
}
}


\author{Stephane CHAMPELY}

\seealso{plotCor, plotBA, plotMcNeil}

\examples{
data(PrisonStress)
plotSliding(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter")
plotSliding(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter",xlab="Before",ylab="After")

}
