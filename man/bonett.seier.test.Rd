\name{bonett.seier.test}
\alias{bonett.seier.test}
\title{Bonett-Seier test of scale for paired samples}
\description{
Robust test of scale for paired samples based on the mean absolute deviations.
}
\usage{
bonett.seier.test(x, y = NULL, alternative = c("two.sided", "less", "greater"),
         omega = 1, conf.level = 0.95)
}
\arguments{
  \item{x}{First sample}
  \item{y}{Second sample}
  \item{alternative}{Alternative hypothesis}
  \item{omega}{A priori ratio of means absolute deviations}
  \item{conf.level}{Confidence level}
}
\value{
A list with class "htest" containing the following components: 
\item{statistic}{The value of the t-statistic} 
\item{p.value}{The p-value for the test} 
\item{conf.int}{A confidence interval for the ratio of means absolute deviations appropriate to the specified alternative hypothesis} 
\item{estimate}{The estimated means absolute deviations} 
\item{null.value}{The specified hypothesized value of the ratio of means absolute deviations} 
\item{alternative}{A character string describing the alternative hypothesis}
\item{method}{A character string indicating what type of test was performed} 
\item{data.name}{A character string giving the name(s) of the data}

}
\references{
D.G. Bonett and E. Seier. Statistical inference for a ratio of dispersions using paired samples. Journal of Educational
and Behavioral Statistics, 28, 21-30, 2003.
}
\author{Stephane CHAMPELY}


\seealso{var.test, pitman.morgan.test, grambsch.test}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-(rnorm(20)+z)*2
bonett.seier.test(x,y)
}
\keyword{htest}
