\name{Iron}

\docType{data}

\alias{Iron}

\title{Iron data from Preece (1982, Table 10)}

\description{
This dataset presents 10 paired data corresponding to percentages of iron found in compounds with the help of two different methods (take a guess: A & B). It is quite intersting to study rounding effect on hypothesis test (have a look at the examples section).
}

\usage{data(Iron)}

\format{A dataframe with 10 rows and 3 columns: 
\tabular{rlll}{
[,1] \tab Compound \tab factor \tab \cr
[,2] \tab Method_A \tab numeric \tab percentage of iron \cr
[,3] \tab Method_B  \tab numeric \tab percentage of iron \cr
}
}

\source{
Chatfield, C. (1978) Statistics for Technology: A Course in Applied Statistics, 2nd ed. Chapman and Hall: London.
}

\references{
Preece, D.A. (1982) t is for trouble (and textbooks): a critique of some examples of the paired-samples t-test. The Statistician, 31 (2), 169-195.
}

\examples{
data(Iron)

# Visualizing, very nice correlation
# Is this an agreement problem or a comparison problem?
with(Iron,plot(paired(Method_A,MethodB)))

# Significant... p=0.045
with(Iron,t.test(paired(Method_A,MethodB)))

# Looking at data, rounded at 0.1 so they can be +0.05 or -0.05
show(Iron)

# Thus the differences can be +0.1 or -0.1
# Influence of rounding on the t-statistic
with(Iron,t.test(Method_A-MethodB+0.1))
with(Iron,t.test(Method_A-MethodB-0.1))
}


\keyword{datasets}
