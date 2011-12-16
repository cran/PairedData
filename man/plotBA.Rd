\name{plotBA}
\alias{plotBA}
\title{Bland-Altman plot}
\description{
Produce a Bland-Altman plot for paired data, including a confidence region for the mean of the differences.
}
\usage{
plotBA(df,condition1,condition2,groups=NULL,facet=TRUE,...)
}

\arguments{
  \item{df}{Data frame}
  \item{condition1}{name of the variable corresponding to measurements in the "first" condition}
 \item{condition2}{name of the variable corresponding to measurements in the "second" condition}
 \item{groups}{name of the variable corresponding to groups}
\item{facet}{In case of groups, faceting (TRUE) or overplotting (FALSE)}
 \item{...}{Arguments to be passed to methods, such as graphical parameters (see par)}
}


\references{
\itemize{
\item J.M. Bland and D.G. Altman. Measuring agreement in method comparison studies. Statistical Methods in Medical Research, 8, 135-160, 1999.
\item D.M. Meek. Two macros for producing graphs to assess agreement between two variables. In Proceedings of Midwest SAS Users Group
Annual Meeting, October 2007.
}
}

\author{Stephane CHAMPELY}

\seealso{tmd}

\examples{
data(PrisonStress)
plotBA(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter")
plotBA(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter",groups="Group")
plotBA(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter",groups="Group",facet=FALSE)
}
