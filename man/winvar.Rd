\name{winvar}
\alias{winvar}
\title{Winsorized variance}
\description{
Compute the Winsorized variance for the data in the vector x. tr is the amount of Winsorization which defaults to .2.
(from Wilcox, 2005)
}
\usage{
winvar(x,tr=.2,na.rm=FALSE)
}
\arguments{
  \item{x}{Sample data}
  \item{tr}{Percentage of trimming}
  \item{na.rm}{Missing values}
}
\value{
The winsorized variance
}
\references{Wilcox, R.R. (2005). Introduction to robust estimation and hypothesis testing. Academic Press.}
\author{Stephane CHAMPELY (from Wilcox, 2005)}


\seealso{var}

\examples{
x<-rnorm(20)
winvar(x)
sqrt(winvar(x))/0.6421
}
