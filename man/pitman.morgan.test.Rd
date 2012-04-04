\name{pitman.morgan.test}
\alias{pitman.morgan.test}
\alias{pitman.morgan.test.default}
\alias{pitman.morgan.test.paired}

\title{
Pitman-Morgan test of variance for paired samples}

\description{
Classical test of variance for paired samples based on the bivariate normal distribution.
}
\usage{
pitman.morgan.test(x, ...)

\method{pitman.morgan.test}{default}(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
    omega = 1, conf.level = 0.95,...)

\method{pitman.morgan.test}{paired}(x, ...)
}

\arguments{
  \item{x}{
first sample or an object of class paired.}
  \item{y}{
second sample.}
  
  \item{alternative}{
alternative hypothesis.}
  \item{omega}{
a priori ratio of variances.}

  \item{conf.level}{
confidence level.}

  \item{\dots}{
further arguments to be passed to or from methods.}
}
\value{
A list with class "htest" containing the following components: 
\item{statistic}{the value of the F-statistic.} 
\item{parameter}{the degrees of freedom for the F-statistic.} 
\item{p.value}{the p-value for the test.} 
\item{conf.int}{a confidence interval for the ratio of variances appropriate to the specified alternative hypothesis.} 
\item{estimate}{the estimated variances.} 
\item{null.value}{the specified hypothesized value of the ratio of variances.} 
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating what type of test was performed.} 
\item{data.name}{a character string giving the name(s) of the data.}

}
\references{
\itemize{
\item W. A. Morgan. A test for the significance of the difference between two variances in a sample from a normal bivariate distribution. Biometrika, 31, 13-19, 1939.
\item E.J.G. Pitman. A note on normal correlation. Biometrika, 31, 9-12, 1939.
}
}

\author{Stephane CHAMPELY}


\seealso{var.test, bonett.seier.test, grambsch.test}

\examples{
data(HORSEbegin)
pitman.morgan.test(HORSEbegin$Actual,HORSEbegin$Imaginary)

p<-with(HORSEbegin,paired(Actual,Imaginary))
pitman.morgan.test(p)
}
\keyword{htest}

