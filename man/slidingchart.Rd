\name{slidingchart}
\alias{slidingchart}
\alias{slidingchart,paired-method}

\title{Sliding square plot}

\description{
Draw a sliding square plot for paired data which mixes the usual scatterplot with the tukey mean-difference plot.
}
\usage{
\S4method{slidingchart}{paired}(object,...)
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

\references{
\itemize{
\item Rosenbaum, P.R. (1989) Exploratory plot for paired data. American Statistician, 43, 108-110.
\item Pontius, J.S. and Schantz, R.M. (1994) Graphical analyses of a twoperiod crossover design. The American Statistician, 48, 249-253.
\item Pruzek, R.M. and Helmreich, J.E. (2009) Enhancing dependent sample analyses with graphics. Journal of Statistics Education, 17.
}
}


\author{Stephane CHAMPELY}

\seealso{plot}

\examples{
data(PrisonStress)
with(PrisonStress,slidingchart(paired(PSSbefore,PSSafter)))
}
