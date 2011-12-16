\name{grambsch.test}
\alias{grambsch.test}
\title{Grambsch test of scale for paired samples}
\description{
Robust test of scale for paired samples
}
\usage{
grambsch.test(x, y = NULL, alternative = c("two.sided", "less", "greater"))
}
\arguments{
  \item{x}{First sample}
  \item{y}{Second sample}
  \item{alternative}{Alternative hypothesis}
}
\value{
A list with class "htest" containing the following components: 
\item{statistic}{The value of the F-statistic} 
\item{p.value}{The p-value for the test} 
\item{null.value}{The specified hypothesized value of the ratio of variances (=1!)} 
\item{alternative}{A character string describing the alternative hypothesis}
\item{method}{A character string indicating what type of test was performed} 
\item{data.name}{A character string giving the name(s) of the data}

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
}
\keyword{htest}
