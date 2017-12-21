\name{paired.plotMcNeil}
\alias{paired.plotMcNeil}

\title{Parallel lines plot}
\description{
Produce a parallel lines plot for paired data.
}
\usage{
paired.plotMcNeil(df, condition1, condition2, groups = NULL, subjects,facet = TRUE, ...) 
}

\arguments{
  \item{df}{
a data frame.
}
  \item{condition1}{
name of the variable corresponding to the second sample.
}
  \item{condition2}{
name of the variable corresponding to the first sample.
}

  \item{groups}{
names of the variable corresponding to groups (optional).
}
  \item{subjects}{
names of the variable corresponding to subjects.
}
  \item{facet}{
faceting or grouping strategy for plotting?
}

\item{\dots}{
further arguments to be passed to methods.}
}

\value{
a graphical object of class ggplot.
}

\references{
McNeil, D.R. (1992) On graphing paired data. The American Statistician, 46 :307-310.
}

\author{Stephane CHAMPELY}

\seealso{plotBA}


\examples{
data(PrisonStress)
paired.plotMcNeil(PrisonStress,"PSSbefore","PSSafter",subjects="Subject")
}

