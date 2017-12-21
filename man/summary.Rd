\name{summary}
\alias{summary}
\alias{summary.paired}
\alias{summary,paired-method}

\title{Summary statistics for paired samples}
\description{
Classical and robust statistics (location, scale and correlation) for paired samples.
}
\usage{
\S4method{summary}{paired}(object,tr=0.2)
}

\arguments{
  \item{object}{an object of class paired.} 
  \item{tr}{percenatge of trimming.} 
}

\value{
A list with a first table corresponding to location and scale statistics and a second table to Pearson and winsorized correlation.

The first table contains four rows corresponding to calculations for x, y, x-y and (x+y)/2 variables.
The location and scale statistics are given in columns.
  
\item{n}{sample size.} 
\item{mean}{mean.} 
\item{median}{median.} 
\item{trim}{trimmed mean (tr=0.2)} 
\item{sd}{standard deviation.} 
\item{IQR}{interquartile range (standardised to be consistent with the sd in the normal case)} 
\item{median ad}{median of absolute deviations (standardised)}
\item{mean ad}{mean of absolute deviations (standardised)} 
\item{sd(w)}{winsorised standard deviation (tr=0.2 and standardised)}
\item{min}{minimum value.}
\item{max}{maximum value.}
}

\author{Stephane CHAMPELY}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-rnorm(20)+z+1
p<-paired(x,y)
summary(p)
}

\keyword{methods}
