\name{grambsch.Var.test}
\alias{grambsch.Var.test}
\alias{grambsch.Var.test.default}
\alias{grambsch.Var.test.paired}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Grambsch test of scale for paired samples
}
\description{
Robust test of scale for paired samples.}

\usage{
grambsch.Var.test(x, ...)

\method{grambsch.Var.test}{default}(x, y = NULL, alternative = c("two.sided", "less", "greater"),...)

\method{grambsch.Var.test}{paired}(x, ...)
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

\details{
Denoting s=x+y and d=x-y, the test proposed by Grambsch (1994, and called by the author 'modified Pitman test') is based on the fact that var(x)-var(y)=cov(x+y,x-y)=cov(s,d). The values z=(s-mean(s))(d-mean(d)) can be tested for null expectation using a classical t test in order to compare the two variances. Note that the p value is computed using the normal distribution.
}

\references{
Grambsch,P.M. (1994) Simple robust tests for scale differences in paired data. Biometrika, 81, 359-372.
}
\author{Stephane CHAMPELY}


\seealso{Var.test, bonettseier.Var.test}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-(rnorm(20)+z)*2
grambsch.Var.test(x,y)

p<-paired(x,y)
grambsch.Var.test(p)
}
\keyword{htest}

