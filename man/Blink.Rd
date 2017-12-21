\name{Blink}

\docType{data}

\alias{Blink}

\title{Blink data from Preece (1982, Table 2)}

\description{
This dataset presents paired data corresponding to average blink-rate per minute of 12 subjects in an experiment of a visual motor task. They had to steer a pencil along a moving track. Each subject was tested under two conditions : a straight track and an oscillating one.
Note that the values from subjects 1 and 2 are somewhat different. 
}

\usage{data(Blink)}

\format{A dataframe with 12 rows and 3 columns:
\tabular{rlll}{
[,1] \tab Subject \tab factor \tab \cr
[,2] \tab Straight \tab numeric \tab blink rate in first condition \cr
[,3] \tab Oscillating \tab numeric \tab blink rate in second condition \cr
}
}

\source{
Wetherhill, G.B. (1972) Elementary Statistical Methods, 2nd ed Chapman and Hall: London.
}

\references{
Preece D.A. (1982) t is for trouble (and textbooks): a critique of some examples of the paired-samples t-test. The Statistician, 31 (2), 169-195.
}

\examples{
data(Blink)

# Visualizing two "outliers"
with(Blink,plot(paired(Straight,Oscillating),type="profile"))

# Interestingly, the differences for the two outliers are quite "normal"
# so their influence on the t test is negligible
with(Blink,qqnorm(Straight-Oscillating))
with(Blink,qqline(Straight-Oscillating))
}

\keyword{datasets}
