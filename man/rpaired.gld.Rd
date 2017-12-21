\name{rpaired.gld}
\alias{rpaired.gld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate paired samples
}
\description{
Simulate paired data with a given correlation (Kendall's tau=(2/pi)arcsine(r)) and marginals being
 Generalized Tukey-Lambda (G-TL) distributions.
}

\usage{
rpaired.gld(n, d1=c(0.000,0.1974,0.1349,0.1349), d2=c(0.000,0.1974,0.1349,0.1349), r)
}

\arguments{
  \item{n}{sample size.}
  \item{d1}{vector of four parameters for the first G-TL distribution.}
  \item{d2}{vector of four parameters for the second G-TL distribution.}
\item{r}{correlation.}
}


\value{
An object of class paired.
}


\references{
Grambsch, P.M. (1994) Simple robust tests for scale differences in paired data. Biometrika, 81, 359-372.
}
\author{Stephane CHAMPELY}


\seealso{rpaired.contaminated}

\examples{
rpaired.gld(n=30,r=0.5)

data(lambda.table)
p<-rpaired.gld(n=30,d1=lambda.table[7,],d2=lambda.table[7,],r=0.5)
plot(p)
}