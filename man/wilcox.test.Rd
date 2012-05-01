\name{wilcox.test}
\alias{wilcox.test}
\alias{wilcox.test.paired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wilcoxon's signed rank test for paired data
}
\description{
A method designed for objects of class paired.
}
\usage{
\method{wilcox.test}{paired}(x, ...)
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
\item{statistic}{the value of V statistic.} 
\item{parameter}{the parameter(s) for the exact distribution of the test statistic.} 
\item{p.value}{the p-value for the test.} 
\item{null.value}{the true location shift mu.} 
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating what type of test was performed (always paired here)} 
\item{data.name}{a character string giving the name(s) of the data.}
\item{conf.int}{	a confidence interval for the location parameter. (Only present if argument conf.int = TRUE.)}
\item{estimate}{	an estimate of the location parameter. (Only present if argument conf.int = TRUE.)}
}


\author{
Stephane Champely
}


\seealso{
yuen.test
}

\examples{
data(PrisonStress)
with(PrisonStress,wilcox.test(PSSbefore,PSSafter))
with(PrisonStress,wilcox.test(PSSbefore,PSSafter,paired=TRUE))
with(PrisonStress,wilcox.test(paired(PSSbefore,PSSafter)))
}


