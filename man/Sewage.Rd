\name{Sewage}

\docType{data}

\alias{Sewage}

\title{Chlorinating sewage data from Preece (1982, Table 9)}

\description{
This dataset presents 8 paired data corresponding to log coliform densities per ml for 2 sewage chlorination methods on each of 8 days.}

\usage{data(Sewage)}

\format{A dataframe with 8 rows and 3 columns: 
\tabular{rlll}{
[,1] \tab Day \tab numeric \tab \cr
[,2] \tab Method_A \tab numeric \tab log density \cr
[,3] \tab Method_B \tab numeric \tab log density \cr
}
}

\source{
Wetherill, G.B. (1972) Elementary Statistical Methods, 2nd ed. Chapman and Hall: London.
}

\references{
Preece, D.A. (1982) t is for trouble (and textbooks): a critique of some examples of the paired-samples t-test. The Statistician, 31 (2), 169-195.
}

\examples{
data(Sewage)

# Visualising
with(Sewage,plot(paired(Method_A,Method_B),type="profile"))

# Basic paired t-test
with(Sewage,t.test(paired(Method_A,Method_B)))

# Influence of the 0.1 rounding on the t-test
 with(Sewage,t.test(Method_A-Method_B-0.1))
 with(Sewage,t.test(Method_A-Method_B+0.1))

}


\keyword{datasets}
