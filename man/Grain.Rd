\name{Grain}

\docType{data}

\alias{Grain}

\title{Grain data from Preece (1982, Table 5)}

\description{
This dataset presents 9 paired data corresponding to the grain yields of Great Northern and Big Four oats grown in "adjacent" plots.
}

\usage{data(Grain)}

\format{A dataframe with 9 rows and 3 columns: 
\tabular{rlll}{
[,1] \tab Year \tab factor \tab \cr
[,2] \tab GreatNorthern \tab numeric \tab grain yield (bushels per acre) \cr
[,3] \tab BigFour \tab numeric \tab grain yield \cr
}
}

\source{
LeClerg, E.L., Leonard, W.H. & Clark, A.G. (1962) Field Plot Technique. Burgess: Minneapolis.
}

\references{
Preece, D.A. (1982) t is for trouble (and textbooks): a critique of some examples of the paired-samples t-test. The Statistician, 31 (2), 169-195.
}

\examples{
data(Grain)

# Usual visualization for paired data (2 clusters?)
with(Grain, plot(paired(GreatNorthern,BigFour)))

# Are they actually "adjacent" plots? 
# Why this variable Year?
# Is there any time trend?
with(Grain, plot(Year,GreatNorthern,type="o"))
with(Grain, plot(Year,BigFour,type="o"))

}


\keyword{datasets}
