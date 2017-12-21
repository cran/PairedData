\name{Meat}

\docType{data}

\alias{Meat}

\title{Meat data from Preece (1982, Table 4)}

\description{
This dataset presents 20 paired data corresponding to the percentage of fat in samples of meat using two different methods: AOAC and Babcock.
}

\usage{data(Meat)}

\format{A dataframe with 20 rows and 3 columns: 
\tabular{rlll}{
[,1] \tab AOAC \tab numeric \tab percentage of fat \cr
[,2] \tab Babcock \tab numeric \tab percentage of fat \cr
[,3] \tab MeatType \tab factor \tab meat type \cr
}
}

\source{
Tippett, L.H.C. (1952) Technological Applications of Statistics. Williams and Norgate: London.
}

\references{
Preece, D.A. (1982) t is for trouble (and textbooks): a critique of some examples of the paired-samples t-test. The Statistician, 31 (2), 169-195.
}

\examples{
data(Meat)

# Presence of clusters or...
with(Meat,plot(paired(AOAC,Babcock)))

# group effect according to Meat type?
with(Meat,plot(paired(AOAC,Babcock),group=MeatType))
with(Meat,plot(paired(AOAC,Babcock),group=MeatType,facet=FALSE))
}

\keyword{datasets}
