\name{rpaired.contaminated}
\alias{rpaired.contaminated}
\title{Simulate paired samples}
\description{
Simulate paired data with a given correlation (Kendall's tau=(2/pi)arcsine(r)) and 
marginals being contaminated normal distributions: (1-e)F(x)+eF(x/K) where F is the cumulative standard normal distribution, e the percentage of contamination and
K a scale parameter.
}
\usage{
rpaired.contaminated(n,d1=c(0.1,10,1),d2=c(0.1,10,1),r=0.5)
}
\arguments{
  \item{n}{Sample size}
  \item{d1}{Vector of 3 parameters for the first contaminated normal distribution}
  \item{d2}{Vector of 3 parameters for the second contaminated normal distribution}
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

\seealso{rpaired.gld}

\examples{
rpaired.contaminated(n=30,r=0.25)
}
