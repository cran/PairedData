\name{rpaired.gld}
\alias{rpaired.gld}
\title{Simulate paired samples}
\description{
Simulate paired data with a given correlation (Kendall's tau=(2/pi)arcsine(r)) and marginals being
 Generalized Tukey-Lambda distributions.
}
\usage{
rpaired.gld(n,d1,d2,r)
}
\arguments{
  \item{n}{Sample size}
  \item{d1}{Vector of 4 parameters for the first G-TL distribution}
  \item{d2}{Vector of 4 parameters for the second G-TL distribution}
\item{r}{Correlation}
}
\value{
A list with two components: 
\item{x}{First sample} 
\item{y}{Second sample} 
}
\references{
P.M. Grambsch. Simple robust tests for scale differences in paired data. Biometrika, 81, 359-372, 1994.
}
\author{Stephane CHAMPELY}


\seealso{rpaired.contaminated}

\examples{
rpaired.gld(n=30,d1=c(0,0.1974,0.1349,0.1349),d2=c(0,0.1974,0.1349,0.1349),r=0.5)
}
