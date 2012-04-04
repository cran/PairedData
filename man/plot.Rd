\name{plot}
\alias{plot.paired}
\alias{plot,paired-method}



\title{~~ Methods for Function \code{plot} ~~ }
\description{
Plot an object of class paired.
}
\usage{
\S4method{plot}{paired}(x, groups=NULL, facet=TRUE,...)
}

\arguments{

  \item{x}{a paired object created by the
    \code{paired} function.}
  \item{groups}{a factor (optional).}
  \item{facet}{faceting or grouping strategy for plotting?.}
  \item{\dots}{arguments to be passed to methods.}

}

\value{
an graphical object of class ggplot.
}

\examples{
data(HORSEbegin)
p<-with(HORSEbegin,paired(Actual,Imaginary))
plot(p)
}

\keyword{methods}
