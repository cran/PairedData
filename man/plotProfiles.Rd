\name{plotProfiles}
\alias{plotProfiles}
\alias{plotProfiles,paired-method}

\title{Profile plot}
\description{
Produce a profile plot or before-after plot  or 1-1 plot for paired data.
}
\usage{
\S4method{plotProfiles}{paired}(object,groups=NULL,subjects=NULL,facet=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
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
arguments to be passed to methods.
}
}


\value{
A graphical object of class ggplot.
}


\references{N.J. Cox. Speaking data: graphing agreement and disagreement. The Stata Journal, 4, 329-349, 2004.}


\author{Stephane CHAMPELY}

\seealso{plotBA,plotMcNeil}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-rnorm(20)+z+1
p<-paired(x,y)
plotProfiles(p)
}