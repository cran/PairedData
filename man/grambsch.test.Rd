\name{grambsch.test}
\alias{grambsch.test}
\alias{grambsch.test.default}
\alias{grambsch.test.paired}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Grambsch test of scale for paired samples
}
\description{
Robust test of scale for paired samples.}

\usage{
grambsch.test(x, ...)

\method{grambsch.test}{default}(x, y = NULL, alternative = c("two.sided", "less", "greater"),...)

\method{grambsch.test}{paired}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
first sample or an object of class paired.
}
  \item{y}{
second sample.
}
  
  \item{alternative}{
alternative hypothesis.}
  \item{\dots}{
further arguments to be passed to or from methods.}
}
\value{
A list with class "htest" containing the following components: 
\item{statistic}{the value of the F-statistic.} 
\item{p.value}{the p-value for the test.} 
\item{null.value}{the specified hypothesized value of the ratio of variances (=1!)} 
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating what type of test was performed.} 
\item{data.name}{a character string giving the name(s) of the data.}

}
\references{
P.M. Grambsch. Simple robust tests for scale differences in paired data. Biometrika, 81, 359-372, 1994.
}
\author{Stephane CHAMPELY}


\seealso{var.test, pitman.morgan.test, bonett.seier.test}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-(rnorm(20)+z)*2
grambsch.test(x,y)

p<-paired(x,y)
grambsch.test(p)
}
\keyword{htest}

