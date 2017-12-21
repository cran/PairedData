\name{Grain2}

\docType{data}

\alias{Grain2}

\title{Wheat grain data from Preece (1982, Table 12)}

\description{
This dataset presents 6 paired data corresponding to the grain yields of two wheat varieties grown on pairs of plots.}

\usage{data(Grain2)}

\format{A dataframe with 6 rows and 3 columns: 
\tabular{rlll}{
[,1] \tab Plot \tab factor \tab \cr
[,2] \tab Variety_1  \tab numeric \tab grain yield (bushels per acre) \cr
[,3] \tab Variety_2 \tab numeric \tab grain yield \cr
}
}

\source{
Balaam, N.L. (1972) Fundamentals of Biometry. The Science of Biology Series (ed J.D. Carthy and J.F. Sutcliffe), No3, Allen and Unwin: London.
}

\references{
Preece, D.A. (1982) t is for trouble (and textbooks): a critique of some examples of the paired-samples t-test. The Statistician, 31 (2), 169-195.
}

\examples{
data(Grain2)

# A very small data set
print(Grain2)

# The paired t test is the test of the differences
with(Grain2,t.test(Variety_1,Variety_2,paired=TRUE))
with(Grain2,t.test(Variety_1-Variety_2))

# The data are actually rounded to the nearest integer
# So they can be somewhere between +0.5 or -0.5
# and thus the differences between +1 or -1
# The possible t values can be simulated by:
simulating.t<-numeric(1000)
for(i in 1:1000){
simulating.t[i]<-with(Grain2,t.test(Variety_1-Variety_2+runif(6,-1,1)))$stat
}
hist(simulating.t)
abline(v=with(Grain2,t.test(Variety_1-Variety_2))$stat,lty=2)
}

\keyword{datasets}
