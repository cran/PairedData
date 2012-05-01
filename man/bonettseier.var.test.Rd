\name{bonettseier.var.test}
\alias{bonettseier.var.test}
\alias{bonettseier.var.test.default}
\alias{bonettseier.var.test.paired}

\title{
Bonett-Seier test of scale for paired samples
}
\description{
Robust test of scale for paired samples based on the mean absolute deviations.
}
\usage{
bonettseier.var.test(x, ...)

\method{bonettseier.var.test}{default}(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
omega = 1, conf.level = 0.95,...)

\method{bonettseier.var.test}{paired}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
first sample or object of class paired.
}
  \item{y}{
second sample.}
  
  \item{alternative}{
alternative hypothesis.}
  \item{omega}{
a priori ratio of means absolute deviations.}

  \item{conf.level}{
confidence level.}

  \item{\dots}{
further arguments to be passed to or from methods.}
}
\value{
A list with class "htest" containing the following components: 
\item{statistic}{the value of the z-statistic.} 
\item{p.value}{the p-value for the test.} 
\item{conf.int}{a confidence interval for the ratio of means absolute deviations appropriate to the specified alternative hypothesis.} 
\item{estimate}{the estimated means absolute deviations.} 
\item{null.value}{the specified hypothesized value of the ratio of means absolute deviations.} 
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating what type of test was performed.} 
\item{data.name}{a character string giving the name(s) of the data.}

}

\references{
Bonett, D.G. and Seier E. (2003) Statistical inference for a ratio of dispersions using paired samples. Journal of Educational
and Behavioral Statistics, 28, 21-30.
}

\author{Stephane CHAMPELY}


\seealso{var.test, grambsch.var.test}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-(rnorm(20)+z)*2
bonettseier.var.test(x,y)

data(anscombe2)
p<-with(anscombe2,paired(X1,Y1))
bonettseier.var.test(p)
}
\keyword{htest}

