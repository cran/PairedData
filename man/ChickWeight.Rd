\name{ChickWeight}

\docType{data}

\alias{ChickWeight}

\title{Chick weight data from Preece (1982, Table 11)}

\description{
This dataset presents 10 paired data corresponding to the weights of chicks, two from ten families, reared in confinement or on open range.}

\usage{data(ChickWeight)}

\format{A dataframe with 10 rows and 3 columns:
\tabular{rlll}{
[,1] \tab Chicks \tab factor \tab \cr
[,2] \tab Confinement \tab numeric \tab chick weight (ounces) \cr
[,3] \tab OpenRange \tab numeric \tab chick weight \cr
}
}

\source{
Paterson, D.D. (1939) Statistical Techniques in Agricultural Research. McGrw-Hill: New-York.
}

\references{
Preece, D.A. (1982) t is for trouble (and textbooks): a critique of some examples of the paired-samples t-test. The Statistician, 31 (2), 169-195.
}

\examples{
data(ChickWeight)

# Look at the interesting discussion in Preece (1982)
# about degree of precision and t test
with(ChickWeight,plot(paired(Confinement,OpenRange)))
with(ChickWeight,stem(Confinement-OpenRange,scale=2))
}


\keyword{datasets}
