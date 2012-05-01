\name{paired.plotBA}
\alias{paired.plotBA}

\title{Bland-Altman plot}
\description{
Produce a Bland-Altman plot for paired data, including a confidence region for the mean of the differences.
}
\usage{
paired.plotBA(df, condition1, condition2, groups = NULL, 
facet = TRUE, ...)
}

\arguments{
  \item{df}{
a data.frame.
}
  \item{condition1}{
name of the variable corresponding to the first sample.
}
  \item{condition2}{
name of the variable corresponding to the first sample.
}
\item{groups}{
name of the variable corresponding to the groups (optional).
}

  \item{facet}{
faceting or grouping strategy for plotting?
}

\item{\dots}{
arguments to be passed to methods}
}

\value{
a graphical object of class ggplot.
}

\references{
\itemize{
\item Bland, J.M. and Altman D.G. (1999) Measuring agreement in method comparison studies. Statistical Methods in Medical Research, 8, 135-160.
\item Meek, D.M. (2007) Two macros for producing graphs to assess agreement between two variables. In Proceedings of Midwest SAS Users Group Annual Meeting, October 2007.
}
}


\author{Stephane CHAMPELY}

\seealso{tmd}

\examples{
data(PrisonStress)
paired.plotBA(PrisonStress,"PSSbefore","PSSafter")

# Extending the resulting ggplot object by faceting
paired.plotBA(PrisonStress,"PSSbefore","PSSafter")+facet_grid(~Group)
}

