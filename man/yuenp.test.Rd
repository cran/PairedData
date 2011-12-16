\name{yuenp.test}
\alias{yuenp.test}
\title{Yuen test for paired samples}
\description{
Robust test of location for paired samples based on trimmed means, greatly inspired from Wilcox (2005)
}

\usage{
yuenp.test(x, y = NULL, alternative = c("two.sided", "less", "greater"),
         mu = 0, conf.level = 0.95,tr=0.2)
}

\arguments{
  \item{x}{First sample}
  \item{y}{Second sample}
  \item{alternative}{Alternative hypothesis}
  \item{mu}{A priori difference of trimmed means}
  \item{conf.level}{Confidence level}
  \item{tr}{Percentage of trimming}
}

\value{
A list with class "htest" containing the following components: 
\item{statistic}{The value of the t-statistic.} 
\item{parameter}{The degrees of freedom for the t-statistic.} 
\item{p.value}{The p-value for the test.} 
\item{conf.int}{A confidence interval for the difference in trimmed means appropriate to the specified alternative hypothesis.} 
\item{estimate}{The estimated trimmed means} 
\item{null.value}{The specified hypothesized value of the trimmed mean difference.} 
\item{alternative}{A character string describing the alternative hypothesis.}
\item{method}{A character string indicating what type of test was performed.} 
\item{data.name}{A character string giving the name(s) of the data.}

}

\references{Wilcox, R.R. (2005). Introduction to robust estimation and hypothesis testing. Academic Press.}

\author{Stephane CHAMPELY}


\seealso{t.test}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-rnorm(20)+z+1
yuenp.test(x,y)
}
\keyword{htest}
