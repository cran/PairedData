\name{plotBA}
\alias{plotBA}
\alias{plotBA,paired-method}

\title{Bland-Altman plot}
\description{
Produce a Bland-Altman plot for paired data, including a confidence region for the mean of the differences.
}
\usage{
\S4method{plotBA}{paired}(object,groups=NULL,facet=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class paired.
}
  \item{groups}{
factor corresponding to groups (optional).
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
\item J.M. Bland and D.G. Altman. Measuring agreement in method comparison studies. Statistical Methods in Medical Research, 8, 135-160, 1999.
\item D.M. Meek. Two macros for producing graphs to assess agreement between two variables. In Proceedings of Midwest SAS Users Group
Annual Meeting, October 2007.
}
}


\author{Stephane CHAMPELY}

\seealso{tmd}

\examples{
data(PrisonStress)
with(PrisonStress,plotBA(paired(PSSbefore,PSSafter)))
with(PrisonStress,plotBA(paired(PSSbefore,PSSafter),
groups=Group))
with(PrisonStress,plotBA(paired(PSSbefore,PSSafter),
groups=Group,facet=FALSE))
}

