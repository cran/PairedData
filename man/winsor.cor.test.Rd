\name{winsor.cor.test}
\alias{winsor.cor.test}
\alias{winsor.cor.test.default}
\alias{winsor.cor.test.paired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Winsorized correlation test (for paired data)
}
\description{
Test for association between paired samples, using winsorized correlation coefficient.
}
\usage{
winsor.cor.test(x, ...)

\method{winsor.cor.test}{default}(x, y, tr=0.2,alternative = c("two.sided", "less", "greater"), ...)


\method{winsor.cor.test}{paired}(x,tr=0.2,alternative = c("two.sided", "less", "greater"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class paired or the first variable.
}
  \item{y}{
second variable.
}

\item{tr}{
percentage of winsorizing.
}
\item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
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
\item{estimate}{the winsorized correlation.} 
\item{null.value}{the specified hypothesized value of the winsorized correlation (=0).} 
\item{alternative}{a character string describing the alternative hypothesis.}
\item{data.name}{a character string giving the name(s) of the data.}
}

\author{
Stephane Champely
}

\seealso{cor.test}

\examples{
data(PrisonStress)
with(PrisonStress,winsor.cor.test(PSSbefore,PSSafter))
with(PrisonStress,winsor.cor.test(paired(PSSbefore,PSSafter)))
}


