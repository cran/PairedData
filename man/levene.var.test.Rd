\name{levene.var.test}
\alias{levene.var.test}
\alias{levene.var.test.default}
\alias{levene.var.test.paired}

\title{
Levene test of scale for paired samples
}
\description{
Robust test of scale for paired samples based on absolute deviations from the trimmed means (or medians), called extended Brown-Forsythe test in Wilcox (1989).
}
\usage{
levene.var.test(x, ...)

\method{levene.var.test}{default}(x, y = NULL,
       alternative = c("two.sided", "less", "greater"),
       mu = 0,conf.level = 0.95,location=c("trim","median"),
tr=0.1, ...)


\method{levene.var.test}{paired}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
first sample or object of class paired.
}
  \item{y}{
second sample.}
  
  \item{alternative}{
alternative hypothesis.}
  
\item{mu}{
the location parameter mu.}

  \item{conf.level}{
confidence level.}

  \item{location}{
location parameter for centering: trimmed mean or median.}

  \item{tr}{
percentage of trimming.}


  \item{\dots}{
further arguments to be passed to or from methods.}
}

\value{
A list with class "htest" containing the components of a paired t test.
}

\details{
The data are transformed as deviations from the trimmed mean: X=abs(x-mean(x,tr=0.1)) and Y=(y-mean(y,tr=0.1)). A paired t test is then carried out on X and Y.
}

\references{
Wilcox, R.R. (1989) Comparing the variances of dependent groups. Psychometrika, 54, 305-315.
}

\author{Stephane CHAMPELY}


\seealso{var.test, grambsch.var.test}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-(rnorm(20)+z)*2
levene.var.test(x,y)

# Some variations
levene.var.test(x,y,tr=0.2)
levene.var.test(x,y,location="median")


data(anscombe2)
p<-with(anscombe2,paired(X2,Y2))
levene.var.test(p)



}
\keyword{htest}

