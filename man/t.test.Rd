\name{t.test}
\alias{t.test}
\alias{t.test.paired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Student's test test for paired data
}
\description{
A method designed for objects of class paired.
}
\usage{
\method{t.test}{paired}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class paired.
}
  \item{\dots}{
further arguments to be passed to or from methods.
}
}


\value{
A list with class "htest" containing the following components: 
\item{statistic}{the value of the t-statistic.} 
\item{parameter}{the degrees of freedom for the t-statistic.} 
\item{p.value}{the p-value for the test.} 
\item{conf.int}{a confidence interval for the mean appropriate to the specified alternative hypothesis.} 
\item{estimate}{the estimated difference in mean.} 
\item{null.value}{the specified hypothesized value of mean difference.} 
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating what type of test was performed (always paired here)} 
\item{data.name}{a character string giving the name(s) of the data.}
}


\author{
Stephane Champely
}


\seealso{
yuen.test
}

\examples{
data(PrisonStress)
with(PrisonStress,t.test(paired(PSSbefore,PSSafter)))
}


