\name{Corn}

\docType{data}

\alias{Corn}

\title{Corn data (Darwin)}

\description{
This dataset presents 15 paired data corresponding to the final height of corn data (Zea Mays), one produced by cross-fertilization and the other by self-fertilization. These data were used by Fisher (1936) and were published in Andrews and Herzberg (1985).}

\usage{data(Corn)}

\format{A dataframe with 15 rows and 4 columns: 
\tabular{rlll}{
[,1] \tab pair \tab numeric \tab \cr
[,2] \tab pot \tab numeric \tab \cr
[,3] \tab Crossed \tab numeric \tab plant height (inches)\cr
[,4] \tab Self \tab numeric \tab plant height \cr
}
}

\source{
Darwin, C. (1876). The Effect of Cross- and Self-fertilization in the Vegetable Kingdom, 2nd Ed. London: John Murray.
}

\references{
\itemize{
\item Andrews, D. and Herzberg, A. (1985) Data: a collection of problems from many fields for the student and research worker. New York: Springer. 
\item Fisher, R.A. (1936) The design of Experiments. Oliver & Boyd: London}
}

\examples{
data(Corn)

# Visualizing two outliers
with(Corn,slidingchart(paired(Crossed,Self)))

# Very bad matching in these data
with(Corn,cor.test(Crossed,Self))
with(Corn,winsor.cor.test(Crossed,Self))


# So the two-sample test is slightly 
# more interesting than the paired test
with(Corn,t.test(Crossed,Self,var.equal=TRUE))
with(Corn,t.test(Crossed,Self,paired=TRUE))

# The Pitman-Morgan test is influenced by the two outliers
with(Corn,Var.test(paired(Crossed,Self)))
with(Corn,grambsch.Var.test(paired(Crossed,Self)))
with(Corn,bonettseier.Var.test(paired(Crossed,Self)))

# Lastly, is there a pot effect?
with(Corn,plot(paired(Crossed,Self)))
with(Corn,plot(paired(Crossed,Self),group=pot))
}


\keyword{datasets}
