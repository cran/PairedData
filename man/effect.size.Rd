\name{effect.size}
\alias{effect.size}
\alias{effect.size,paired-method}
\title{Effect size computations for paired data}
\description{
Robust and classical effects sizes for paired samples of the form: (Mx-My)/S where Mx and My are location parameters for each sample and S is a scale parameter
}
\usage{
\S4method{effect.size}{paired}(object,tr=0.2)
}

\arguments{
  \item{object}{an object of class paired}
\item{tr}{percentage of trimming}
}

\value{
A table with two rows corresponding to classical (means) and robust (trimmed means, tr=0.2) delta-type effect sizes. The four columns correspond to:

\item{Average}{Numerator is the difference in (trimmed) means, denominator is the average of the two (winsorised and rescaled to be consistent with the standard deviation when the distribution is normal) standard deviations} 
\item{Single (x)}{Denominator is the (winsorised and rescaled) standard deviation of the first sample} 
\item{Single (y)}{Denominator is the (winsorised and rescaled) standard deviation of the second sample} 
\item{Difference}{Numerator is the (trimmed) mean and denominator the (winsorised and rescaled) standard deviation of the differences (x-y)} 
}

\author{Stephane CHAMPELY}

\references{Algina, J., Keselman, H.J. and Penfield, R.D. (2005) Effects sizes and their intervals: 
the two-level repeated measures case. Educational and Psychological Measurement, 65, 241-258.}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-rnorm(20)+z+1
p<-paired(x,y)
effect.size(p)
}
