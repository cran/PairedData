\name{yuen.test}
\alias{yuen.test}
\alias{yuen.test.default}
\alias{yuen.test.formula}
\alias{yuen.test.paired}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Yuen's trimmed mean test
}
\description{
Yuen's test for one, two or paired samples.
}
\usage{
yuen.test(x, ...)

\method{yuen.test}{default}(x, y = NULL, tr = 0.2, alternative = c("two.sided", "less", "greater"),
mu = 0, paired = FALSE, conf.level = 0.95, ...)

\method{yuen.test}{formula}(formula, data, subset, na.action, ...)

\method{yuen.test}{paired}(x, ...)
}

\arguments{
  \item{x}{
first sample or object of class paired.
}
  \item{y}{
second sample.
}
  \item{tr}{
percentage of trimming.
}
  \item{alternative}{
alternative hypothesis.
}
  \item{mu}{
a number indicating the true value of the trimmed mean (or difference in trimmed means if you are performing a two sample test).}
  \item{paired}{
a logical indicating whether you want a paired yuen's test.
}
  \item{conf.level}{
confidence level.
}
  \item{formula}{
a formula of the form y ~ f where y is a numeric variable giving the data values and f a factor with TWO levels giving the corresponding groups.
}
  \item{data}{
an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).
}
  \item{subset}{
an optional vector specifying a subset of observations to be used.
}
  \item{na.action}{
a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").
}

  \item{\dots}{
further arguments to be passed to or from methods.
}
}


\value{
A list with class "htest" containing the following components: 
\item{statistic}{the value of the t-statistic.} 
\item{parameter}{the degrees of freedom for the t-statistic.} 
\item{p.value}{the p-value for the test.} 
\item{conf.int}{a confidence interval for the trimmed mean appropriate to the specified alternative hypothesis.} 
\item{estimate}{the estimated trimmed mean or difference in trimmed means depending on whether it was a one-sample test or a two-sample test.
} 
\item{null.value}{the specified hypothesized value of the trimmed mean or trimmed mean difference depending on whether it was a one-sample test or a two-sample test.} 
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating what type of test was performed.} 
\item{data.name}{a character string giving the name(s) of the data.}

}


\references{
\itemize{
\item Wilcox, R.R. (2005). Introduction to robust estimation and hypothesis testing. Academic Press.
\item Yuen K.K. (1974) The two-sample trimmed t for unequal population variances. Biometrika, 61, 165-170.
}
}

\author{Stephane CHAMPELY}


\seealso{t.test}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-rnorm(20)+z+1

# two-sample test
yuen.test(x,y)

# one-sample test
yuen.test(x)

# paired-sample tests
yuen.test(x,y,paired=TRUE)

p<-paired(x,y)
yuen.test(p)
}

\keyword{htest}