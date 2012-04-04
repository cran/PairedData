\name{plotMcNeil}
\alias{plotMcNeil}
\alias{plotMcNeil,paired-method}

\title{Parallel lines plot}
\description{
Produce a parallel lines plot for paired data.
}
\usage{
\S4method{plotMcNeil}{paired}(object,groups=NULL,subjects=NULL,facet=TRUE,...)
}

\arguments{
  \item{object}{
an object of class paired.
}
  \item{groups}{
factor corresponding to groups (optional).
}
  \item{subjects}{
vector corresponding to subjects (optional).
}
  \item{facet}{
faceting or grouping strategy for plotting?
}

\item{\dots}{
arguments to be passed to methods.}
}

\value{
A graphical object of class ggplot.
}

\references{
D.R. McNeil. On graphing paired data. The American Statistician, 46 :307-310, 1992.
}

\author{Stephane CHAMPELY}

\seealso{plotBA}


\examples{
data(PrisonStress)
with(PrisonStress,plotMcNeil(paired(PSSbefore,PSSafter)))
with(PrisonStress,plotMcNeil(paired(PSSbefore,PSSafter),
groups=Group))
}

