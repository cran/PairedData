\name{sandvikolsson.var.test}
\alias{sandvikolsson.var.test}
\alias{sandvikolsson.var.test.default}
\alias{sandvikolsson.var.test.paired}

\title{
Sandvik-Olsson test of scale for paired samples
}
\description{
Robust test of scale for paired samples based on the absolute deviations from the trimmed means (or medians).
}
\usage{
sandvikolsson.var.test(x, ...)

\method{sandvikolsson.var.test}{default}(x, y = NULL,
            alternative = c("two.sided", "less", "greater"),
            mu = 0, exact = NULL, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95,location=c("trim","median"),tr=0.1, ...)

\method{sandvikolsson.var.test}{paired}(x, ...)
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

\item{exact}{	
a logical indicating whether an exact p-value should be computed.
}

\item{correct}{
a logical indicating whether to apply continuity correction in the normal approximation for the p-value.
}

\item{conf.int}{
a logical indicating whether a confidence interval should be computed.
}

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
A list with class "htest" containing the components of a wilcoxon signed-rank test.
}

\details{
The data are transformed as deviations from the trimmed mean: X=abs(x-mean(x,tr=0.1)) and Y=(y-mean(y,tr=0.1)). A wilcoxon signed-rank test is then carried out on X and Y.
}

\references{
Sandvik, L. and Olsson, B. (1982) A nearly distribution-free test for comparing dispersion in paired samples. Biometrika, 69, 484-485.
}

\author{Stephane CHAMPELY}


\seealso{var.test, grambsch.var.test}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-(rnorm(20)+z)*2
sandvikolsson.var.test(x,y)

p<-paired(x,y)
sandvikolsson.var.test(p)


# some variations
sandvikolsson.var.test(p,tr=0.2)
sandvikolsson.var.test(p,location="median")

}
\keyword{htest}

