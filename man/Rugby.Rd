\name{Rugby}

\docType{data}

\alias{Rugby}

\title{Agreement study in rugby expert ratings}

\description{
This dataset gives the ratings on a continuous ten-points scale of two experts about 93 actions
 during several rugby union matches.
}

\usage{data(Rugby)}

\format{A dataframe with 93 rows and 3 columns.
\tabular{rlll}{
[,1] \tab EXPERT.1 \tab numeric \tab First expert ratings \cr
[,2] \tab EXPERT.2 \tab numeric \tab Second expert ratings \cr
[,3] \tab Actions \tab factor \tab Subject label\cr
}
}

\source{Private communication. Mickael Campo, CRIS, Lyon 1 University, FRANCE.}

\examples{
data(Rugby)
with(Rugby,plot(paired(EXPERT.1,EXPERT.2)))
with(Rugby,plot(paired(EXPERT.1,EXPERT.2),type="BA"))
}

\keyword{datasets}
