\name{winval}
\alias{winval}
\title{Data winsorization}
\description{
Compute the Winsorized values for data x. tr is the amount of Winsorization which defaults to .2.
(from Wilcox, 2005)
}
\usage{
winval(x,tr=.2)
}
\arguments{
  \item{x}{Sample data}
  \item{tr}{Percentage of trimming}
}
\value{
The winsorized values
}
\references{Wilcox, R.R. (2005). Introduction to robust estimation and hypothesis testing. Academic Press.}
\author{Stephane CHAMPELY (from Wilcox, 2005)}


\examples{
x<-rnorm(20)
sort(x)
winval(x)
}
