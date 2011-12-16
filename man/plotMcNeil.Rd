\name{plotMcNeil}
\alias{plotMcNeil}
\title{parallel lines plot}
\description{
Produce a parallel lines plot for paired data.
}
\usage{
plotMcNeil(df,condition1,condition2,subjects,groups=NULL,facet=TRUE,...)
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


\references{
D.R. McNeil. On graphing paired data. The American Statistician, 46 :307-310, 1992.
}

\author{Stephane CHAMPELY}

\seealso{plotCor, plotBA}


\examples{
data(PrisonStress)
plotMcNeil(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter",subjects="Subject")
plotMcNeil(df=PrisonStress,condition1="PSSbefore",condition2="PSSafter",subjects="Subject",groups="Group")
}
