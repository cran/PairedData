\name{IceSkating}

\docType{data}

\alias{IceSkating}

\title{Ice skating speed study}

\description{
This dataset gives the speed measurement (m/sec) for seven iceskating dancers using the return leg in flexion
or in extension.
}

\usage{data(IceSkating)}

\format{A dataframe with 7 rows and 3 columns.
\tabular{rlll}{
[,1] \tab Subject \tab factor \tab anonymous subjects \cr
[,2] \tab Extension \tab numeric \tab speed when return leg in extension (m/sec) \cr
[,3] \tab Flexion \tab numeric \tab speed when return leg in flexion (m/sec) \cr
}
}

\source{Private communication. Karine Monteil, CRIS, Lyon 1 University, FRANCE.}

\references{
Haguenauer, M., Legreneur, P., Colloud, F. and Monteil, K.M. (2002) Characterisation of
the Push-off in Ice Dancing: Influence of the Support Leg extension on Performance.
Journal of Human Movement Studies, 43, 197-210.
}

\examples{
data(IceSkating)

# Nothing particular in the paired plot
with(IceSkating,plot(paired(Extension,Flexion),type="McNeil"))

# The differences are normally distributed
with(IceSkating,qqnorm(Extension-Flexion))
with(IceSkating,qqline(Extension-Flexion))

# Usual t test
with(IceSkating,t.test(paired(Extension,Flexion)))
}

\keyword{datasets}
