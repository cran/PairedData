\name{plotSliding}
\alias{plotSliding}
\alias{plotSliding,paired-method}

\title{Sliding square plot}

\description{
Produce a sliding square plot for paired data which mixes the usual scatterplot 
with the tukey mean-difference plot.
}
\usage{
\S4method{plotSliding}{paired}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class paired.
}
\item{\dots}{
arguments to be passed to methods.
}
}


\value{
A graphical object of class ggplot.
}


\references{
\itemize{
\item P.R. Rosenbaum. Exploratory plot for paired data. American Statistician, 43, 108-110, 1989.
\item J.S. Pontius and R.M. Schantz. Graphical analyses of a twoperiod crossover design. The American Statistician, 48, 249-253, 1994.
\item R.M. Pruzek and J.E. Helmreich. Enhancing dependent sample analyses with graphics. Journal of Statistics Education, 17, 2009.
}
}


\author{Stephane CHAMPELY}

\seealso{plotBA,plotProfiles,plotMcNeil}

\examples{
data(PrisonStress)
with(PrisonStress,plotSliding(paired(PSSbefore,PSSafter)))
}
