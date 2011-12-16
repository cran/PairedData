\name{wincor}
\alias{wincor}
\title{Winsorized correlation}
\description{
Compute the Winsorized correlation for paired data x and y. 
tr is the amount of Winsorization which defaults to .2.
(from Wilcox, 2005)
}
\usage{
wincor(x,y,tr=.2)
}
\arguments{
  \item{x}{First sample}
  \item{y}{Second sample}
  \item{tr}{Percentage of trimming}
}
\value{
The winsorized correlation
}
\references{Wilcox, R.R. (2005). Introduction to robust estimation and hypothesis testing. Academic Press.}

\author{Stephane CHAMPELY (from Wilcox, 2005)}


\seealso{cor}

\examples{
z<-rnorm(20)
x<-rnorm(20)+z
y<-rnorm(20)+z
wincor(x,y)
}
