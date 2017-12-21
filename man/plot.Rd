\name{plot}
\alias{plot.paired}
\alias{plot,paired-method}



\title{~~ Methods for Function \code{plot} ~~ }
\description{
Plot an object of class paired.
}
\usage{
\S4method{plot}{paired}(x, groups=NULL,subjects=NULL,

 facet=TRUE,type=c("correlation","BA","McNeil","profile"),...)
}

\arguments{

  \item{x}{a paired object created by the
    \code{paired} function.}
  \item{groups}{a factor (optional).}
  \item{subjects}{subjects name.}
  \item{facet}{faceting or grouping strategy for plotting?}
  \item{type}{type of the plot (correlation, Bland-Altman, McNeil or profile plot).}

  \item{\dots}{arguments to be passed to methods.}

}

\value{
an graphical object of class ggplot.
}

\examples{
data(HorseBeginners)
pd1<-with(HorseBeginners,paired(Actual,Imaginary))
plot(pd1)
plot(pd1,type="BA")
plot(pd1,type="McNeil")
plot(pd1,type="profile")

data(Shoulder)
with(Shoulder,plot(paired(Left,Right),groups=Group))
with(Shoulder,plot(paired(Left,Right),groups=Group,facet=FALSE))
with(Shoulder,plot(paired(Left,Right),
groups=Group,facet=FALSE,type="profile"))+theme_bw()
}

\keyword{methods}
