\name{anscombe2}

\docType{data}

\alias{anscombe2}

\title{Teaching the paired t test}

\description{
This dataset presents four sets of paired samples (n=15), giving the same t statistic (t=2.11) and thus the same
p-value whereas their situations are really diversified (differences in variances, clustering, heteroscedasticity). 
The importance of plotting data is thus stressed.
The name is given from the famous Anscombe's dataset created to study simple linear regression.
}

\usage{data(anscombe2)}

\format{A dataframe with 15 rows, 8 numeric columns of paired data: (X1,Y1) ; (X2,Y2) ; (X3,Y3) ; (X4,Y4), and 1 factor column:  Subjects, giving a label for the subjects.
}

\source{S. Champely, CRIS, Lyon 1 University, FRANCE}

\references{F. Anscombe, Graphs in statistical analysis. The American Statistican, 27, 17-21.}

\examples{
data(anscombe2)
# p=0.05 for the paired t-test
with(anscombe2,plot(paired(X1,Y1),type="BA"))
with(anscombe2,t.test(paired(X1,Y1)))

# Same p but Var(X2)<Var(Y2) and
# correlation in the Bland-Altman plot
with(anscombe2,t.test(paired(X2,Y2)))
with(anscombe2,summary(paired(X2,Y2)))
with(anscombe2,plot(paired(X2,Y2),type="BA"))

# Same p but two clusters
with(anscombe2,plot(paired(X3,Y3),type="BA"))

# Same p but the difference is "linked" to the mean
with(anscombe2,plot(paired(X4,Y4),type="BA"))
}


\keyword{datasets}
