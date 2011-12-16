\name{plotCor}
\alias{plotCor}
\title{Correlation plot}
\description{
Produce a correlation plot for paired data with a 45 degrees reference line.
}
\usage{
plotCor(df,condition1,condition2,groups=NULL,facet=TRUE,...)
}
\arguments{
  \item{df}{Data frame}
  \item{condition1}{name of the variable corresponding to measurements in the "first" condition}
 \item{condition2}{name of the variable corresponding to measurements in the "second" condition}
 \item{groups}{name of the variable corresponding to groups}
\item{facet}{In case of groups, faceting (TRUE) or overplotting (FALSE)}
 \item{...}{Arguments to be passed to methods, such as graphical parameters (see par)}
}

\author{Stephane CHAMPELY}

\seealso{plot}

\examples{
data(PrisonStress)
plotCor(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter")
plotCor(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter",groups="Group")
plotCor(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter",groups="Group",facet=FALSE)
}
