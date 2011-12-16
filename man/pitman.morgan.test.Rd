\name{pitman.morgan.test}
\alias{pitman.morgan.test}
\title{Pitman-Morgan test of variances for paired samples}
\description{
Classical test of variance for paired samples based on the bivariate normal distribution
}
\usage{
pitman.morgan.test(x, y = NULL, alternative = c("two.sided", "less", "greater"),
         omega = 1, conf.level = 0.95)
}
\arguments{
  \item{x}{First sample}
  \item{y}{Second sample}
  \item{alternative}{Alternative hypothesis}
  \item{omega}{A priori ratio of variances}
  \item{conf.level}{Confidence level}
}
\value{
A list with class "htest" containing the following components: 
\item{statistic}{The value of the F-statistic} 
\item{parameter}{The degrees of freedom for the F-statistic} 
\item{p.value}{The p-value for the test} 
\item{conf.int}{A confidence interval for the ratio of variances appropriate to the specified alternative hypothesis} 
\item{estimate}{The estimated variances} 
\item{null.value}{The specified hypothesized value of the ratio of variances} 
\item{alternative}{A character string describing the alternative hypothesis}
\item{method}{A character string indicating what type of test was performed} 
\item{data.name}{A character string giving the name(s) of the data}

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
z<-rnorm(20)
x<-rnorm(20)+z
y<-(rnorm(20)+z)*2
pitman.morgan.test(x,y)
}
\keyword{htest}
