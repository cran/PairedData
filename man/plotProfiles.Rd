\name{plotProfiles}
\alias{plotProfiles}
\title{Profiles plot}
\description{
Produce a profiles plot or before-after plot  or 1-1 plot for paired data.
}
\usage{
plotProfiles(df,condition1,condition2,subjects,groups=NULL,facet=TRUE,...)
}

\arguments{
  \item{df}{Data frame}
  \item{condition1}{name of the variable corresponding to measurements in the "first" condition}
 \item{condition2}{name of the variable corresponding to measurements in the "second" condition}
 \item{subjects}{name of the variable corresponding to subjects}
 \item{groups}{name of the variable corresponding to groups}
\item{facet}{In case of groups, faceting (TRUE) or overplotting (FALSE)}
 \item{...}{Arguments to be passed to methods, such as graphical parameters (see par)}
}

\references{N.J. Cox. Speaking data: graphing agreement and disagreement. The Stata Journal, 4, 329-349, 2004.}


\author{Stephane CHAMPELY}

\seealso{plotCor, plotBA}

\examples{
data(PrisonStress)
plotProfiles(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter",subjects="Subject")
plotProfiles(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter",subjects="Subject")+xlim("PSSbefore","PSSafter")
plotProfiles(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter",subjects="Subject",groups="Group")+xlim("PSSbefore","PSSafter")
plotProfiles(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter",subjects="Subject",groups="Group",facet=FALSE)+xlim("PSSbefore","PSSafter")
}
