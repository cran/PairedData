\name{Barley}

\docType{data}

\alias{Barley}

\title{Barley data from Preece (1982, Table 1)}

\description{
This dataset presents 12 paired data corresponding to the yields of Glabron and Velvet Barley, grown on different farms.
The values from farm 12 are quite different.}

\usage{data(Barley)}

\format{A dataframe with 17 rows and 3 columns:
\tabular{rlll}{
[,1] \tab Farm \tab factor \tab \cr
[,2] \tab Glabron \tab numeric \tab yields (bushels per acre) \cr
[,3] \tab Velvet \tab numeric \tab yields \cr
}

}

\source{
Leonard, W.H. & Clark, A.G. (1939) Field Plot Techniques. Burgess: Minneapolis.
}

\references{
Preece D.A. (1982) t is for trouble (and textbooks): a critique of some examples of the paired-samples t-test. The Statistician, 31 (2), 169-195.
}

\examples{
data(Barley)

# Visualizing a clear outlier
with(Barley,plot(paired(Glabron,Velvet),type="BA"))

# Results form the paired t test and paired Yuen test are similar
with(Barley,t.test(paired(Glabron,Velvet)))
with(Barley,yuen.t.test(paired(Glabron,Velvet)))

# Nevertheless the outlier inflates the location (numerator) and
# scale (denominator) standard statictics for the difference
with(Barley,summary(paired(Glabron,Velvet)))
}

\keyword{datasets}
