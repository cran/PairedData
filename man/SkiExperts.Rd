\name{SkiExperts}

\docType{data}

\alias{SkiExperts}

\title{Actual and imaginary performances in ski}

\description{
This dataset gives the actual and motor imaginary performances (time) in ski for 12 experts.
}

\usage{data(SkiExperts)}

\format{A dataframe with 12 rows and 3 columns.
\tabular{rlll}{
[,1] \tab Subject \tab factor \tab anonymous subjects \cr
[,2] \tab Actual \tab numeric \tab actual performance (sec.)\cr
[,3] \tab Imaginary \tab numeric \tab imaginary performance (sec.)\cr
}
}

\source{Private communication. Aymeric Guillot, CRIS, Lyon 1 University, FRANCE.}

\references{
Louis, M., Collet, C., Champely, S. and Guillot, A. (2012) Differences in motor imagery time when predicting task duration in Alpine skiers and equestrian riders. 
Research Quarterly for Exercise and Sport, 83(1), 86-93.
}

\examples{
data(SkiExperts)

# Visualising
with(SkiExperts,plot(paired(Actual,Imaginary),type="profile"))

# No underestimation of imaginary time for experts
with(SkiExperts,t.test(paired(Actual,Imaginary)))

# But a very interesting increase in dispersion in their
# predicted times
with(SkiExperts,Var.test(paired(Actual,Imaginary)))
}

\keyword{datasets}
