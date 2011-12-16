\name{paired.summary}
\alias{paired.summary}
\title{Summary statistics for paired samples}
\description{
Classical and robust statistics (location and scale) for paired samples}

\usage{
paired.summary(x,y)
}
\arguments{
  \item{x}{First sample}
  \item{y}{Second sample}
}
\value{
A table with four rows corresponding to calculations for x, y, x-y and (x+y)/2 variables.
The location and scale statistics are given in columns.
  
\item{n}{Sample size} 
\item{mean}{Mean} 
\item{median}{Median} 
\item{trim}{Trimmed mean (tr=0.2)} 
\item{sd}{Standard deviation} 
\item{IQR}{Interquartile range (standardised to be consistent with the sd in the normal case)} 
\item{median ad}{Median of absolute deviations (standardised)}
\item{mean ad}{Mean of absolute deviations (standardised)} 
\item{sd(w)}{Winsorised standard deviation (tr=0.2 and standardised)}
\item{min}{Minimum value}
\item{max}{Maximum value}
}

\author{Stephane CHAMPELY}

\seealso{summary}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-rnorm(20)+z+1
paired.summary(x,y)
}
