\name{mcculloch.var.test}
\alias{mcculloch.var.test}
\alias{mcculloch.var.test.default}
\alias{mcculloch.var.test.paired}

\title{
McCulloch test of scale for paired samples
}
\description{
Robust test of scale for paired samples based on spearman coefficient (the default, or kendall or pearson) of the transformed D=x-y and S=x+y.
}
\usage{
mcculloch.var.test(x, ...)

\method{mcculloch.var.test}{default}(x, y = NULL,
       alternative = c("two.sided", "less", "greater"),method= c("spearman","pearson", "kendall"),
       exact = NULL,conf.level = 0.95,continuity = FALSE, ...)

\method{mcculloch.var.test}{paired}(x, ...)
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

\item{method}{
a character string indicating which correlation coefficient is to be used for the test. One of "spearman", "kendall", or "pearson", can be abbreviated.
}
  
\item{exact}{
a logical indicating whether an exact p-value should be computed.}

  \item{conf.level}{
confidence level.}

\item{continuity}{
logical: if true, a continuity correction is used for Spearman's rho when not computed exactly.}

  \item{\dots}{
further arguments to be passed to or from methods.}
}

\value{
A list with class "htest" containing the components of a (Spearman) correlation test.
}


\references{
McCulloch, C.E. (1987) Tests for equality of variances for paired data. Communications in Statistics - Theory and Methods, 16, 1377-1391.}

\author{Stephane CHAMPELY}


\seealso{var.test, grambsch.var.test}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-(rnorm(20)+z)*2
mcculloch.var.test(x,y)

p<-paired(x,y)
mcculloch.var.test(p)

# A variation with kendall tau
mcculloch.var.test(p,method="kendall")

# equivalence with the PitmanMorgan test
mcculloch.var.test(p,method="pearson")
var.test(p)


}
\keyword{htest}

