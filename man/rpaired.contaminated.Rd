\name{rpaired.contaminated}
\alias{rpaired.contaminated}

\title{Simulate paired samples}
\description{
Simulate paired data with a given correlation (Kendall's tau=(2/pi)arcsine(r)) and 
marginals being contaminated normal distributions: (1-eps)*F(x)+eps*F(x/K) where F is the cumulative standard normal distribution, eps the percentage of contamination and
K a scale parameter. Moreover, this marginal can be multiplied by another scale parameter sigma but usually sigma=1.
}

\usage{
rpaired.contaminated(n, d1 = c(0.1, 10, 1), d2 = c(0.1, 10, 1), r = 0.5)
}

\arguments{
  \item{n}{
sample size.
}
  \item{d1}{
vector of 3 parameters for the first contaminated normal distribution (eps,K,sigma).
}
  \item{d2}{
vector of 3 parameters for the second contaminated normal distribution.
}
  \item{r}{
correlation.
}

}



\value{
An object of class paired.
}


\references{
P.M. Grambsch. Simple robust tests for scale differences in paired data. Biometrika, 81, 359-372, 1994.
}

\author{Stephane CHAMPELY}

\seealso{rpaired.gld}

\examples{
rpaired.contaminated(n=30,r=0.25)
}
