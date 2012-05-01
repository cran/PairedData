\name{var.test}
\alias{var.test}
\alias{var.test.default}
\alias{pitman.morgan.test.default}
\alias{var.test.paired}

\title{
Tests of variance(s) for normal distribution(s)}

\description{
Classical tests of variance for one-sample, two-independent samples or paired samples.}
\usage{
\method{var.test}{default}(x, y = NULL, ratio = 1, alternative = c("two.sided", 
    "less", "greater"), paired = FALSE, conf.level = 0.95, ...)

\method{var.test}{paired}(x, ...)

\method{pitman.morgan.test}{default}(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
    ratio = 1, conf.level = 0.95,...)

}

\arguments{
  \item{x}{
first sample or an object of class paired or an object of class lm.}
  \item{y}{
second sample or an object of class lm.}
    \item{ratio}{
a priori ratio of variances (two-samples) or variance (one-sample).}

  \item{alternative}{
alternative hypothesis.}
    \item{paired}{independent (the default) or paired samples.}

  \item{conf.level}{
confidence level.}

  \item{\dots}{
further arguments to be passed to or from methods.}
}
\value{
A list with class "htest" containing the following components: 
\item{statistic}{the value of the X-squared statistic (one-sample) or F-statistic (two-samples).} 
\item{parameter}{the degrees of freedom for the statistic.} 
\item{p.value}{the p-value for the test.} 
\item{conf.int}{a confidence interval for the parameter appropriate to the specified alternative hypothesis.} 
\item{estimate}{the estimated variance(s).} 
\item{null.value}{the specified hypothesized value of the parameter.} 
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating what type of test was performed.} 
\item{data.name}{a character string giving the name(s) of the data.}

}
\references{
\itemize{
\item Morgan, W.A. (1939) A test for the significance of the difference between two variances in a sample from a normal bivariate distribution. Biometrika, 31, 13-19.
\item Pitman, E.J.G. (1939) A note on normal correlation. Biometrika, 31, 9-12.
}
}

\author{Stephane CHAMPELY}


\seealso{bonettseier.var.test, grambsch.var.test}

\examples{
data(HorseBeginners)

#one sample test
var.test(HorseBeginners$Actual,ratio=15)

# two independent samples test
var.test(HorseBeginners$Actual,HorseBeginners$Imaginary)

# two dependent samples test
var.test(HorseBeginners$Actual,HorseBeginners$Imaginary,paired=TRUE)
p<-with(HorseBeginners,paired(Actual,Imaginary))
var.test(p)
}
\keyword{htest}

