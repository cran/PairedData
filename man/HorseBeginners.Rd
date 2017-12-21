\name{HorseBeginners}

\docType{data}

\alias{HorseBeginners}

\title{Actual and imaginary performances in equitation}

\description{
This dataset gives the actual and motor imaginary performances (time) in horse-riding for 8 beginners.
}

\usage{data(HorseBeginners)}

\format{A dataframe with 8 rows and 3 columns.
\tabular{rlll}{
[,1] \tab Subject \tab factor \tab Anonymous subjects \cr
[,2] \tab Actual \tab numeric \tab Actual performance (sec.)\cr
[,3] \tab Imaginary \tab numeric \tab Imaginary performance (sec.)\cr
}
}

\source{Private communication. Aymeric Guillot, CRIS, Lyon 1 University, FRANCE.}

\references{
Louis, M. Collet, C. Champely, S. and Guillot, A. (2010) Differences in motor imagery time when predicting task duration. 
Research Quarterly for Exercise and Sport.
}

\examples{
data(HorseBeginners)

# There is one outlier
with(HorseBeginners,plot(paired(Actual,Imaginary),type="profile"))

# This outlier has a great influence
# on the non robust Pitman-Morgan test of variances
with(HorseBeginners,Var.test(paired(Actual,Imaginary)))
with(HorseBeginners[-1,],Var.test(paired(Actual,Imaginary)))
with(HorseBeginners,grambsch.Var.test(paired(Actual,Imaginary)))
with(HorseBeginners,bonettseier.Var.test(paired(Actual,Imaginary)))
}


\keyword{datasets}
