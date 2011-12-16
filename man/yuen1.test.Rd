\name{yuen1.test}
\alias{yuen1.test}
\title{Yuen test for one sample}
\description{
Robust test of location for one sample based on the trimmed mean, greatly inspired from Wilcox (2005)
}
\usage{
yuen1.test(x,tr=.2,conf.level = 0.95,mu=0,alternative=c("two.sided", "less", "greater"))
}
\arguments{
  \item{x}{Sample}
\item{tr}{Percentage of trimming}
\item{conf.level}{Confidence level}
\item{mu}{A priori value of the (population) trimmed mean}
  \item{alternative}{Alternative hypothesis}
  
  
  
}
\value{
A list with class "htest" containing the following components: 
\item{statistic}{The value of the t-statistic} 
\item{parameter}{The degrees of freedom for the t-statistic} 
\item{p.value}{The p-value for the test} 
\item{conf.int}{A confidence interval for the trimmed mean appropriate to the specified alternative hypothesis} 
\item{estimate}{The estimated trimmed mean} 
\item{null.value}{The specified hypothesized value of the (population) trimmed mean} 
\item{alternative}{A character string describing the alternative hypothesis}
\item{method}{A character string indicating what type of test was performed} 
\item{data.name}{A character string giving the name(s) of the data}

}
\references{Wilcox, R.R. (2005). Introduction to robust estimation and hypothesis testing. Academic Press.}

\author{Stephane CHAMPELY}


\seealso{t.test, yuenp.test, wilcox.test}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-rnorm(20)+z+1
yuen1.test(x-y)
}

\keyword{htest}
