\name{GrapeFruit}

\docType{data}

\alias{GrapeFruit}

\title{Grape Fruit data from Preece (1982, Table 6)}

\description{
This dataset presents paired data corresponding to the percentage of solids recorded in the shaded and exposed halves of 25 grapefruits.}

\usage{data(GrapeFruit)}

\format{A dataframe with 25 rows and 3 columns: 
\tabular{rlll}{
[,1] \tab Fruit \tab numeric \tab \cr
[,2] \tab Shaded \tab numeric \tab percentage of solids in grapefruit \cr
[,3] \tab Exposed \tab numeric \tab percentage of solids \cr
}
}

\source{
Croxton, F.E. & Coxden, D.J. (1955) Applied Genral Statistics, 2nd ed. Chapman and Hall, London.
}

\references{
Preece, D.A. (1982) t is for trouble (and textbooks): a critique of some examples of the paired-samples t-test. The Statistician, 31 (2), 169-195.
}

\examples{
data(GrapeFruit)

# Visualizing a very strange paired distribution
with(GrapeFruit,plot(paired(Shaded,Exposed)))
with(GrapeFruit,plot(paired(Shaded,Exposed),type="BA"))
with(GrapeFruit,plot(paired(Shaded,Exposed),type="McNeil"))
with(GrapeFruit,plot(paired(Shaded,Exposed),type="profile"))

# As underlined by Preece (1982), have a look to
# the distribution of the final digits
show(GrapeFruit)
table(round((GrapeFruit$Shaded*10-floor(GrapeFruit$Shaded*10))*10))
table(round((GrapeFruit$Exposed*10-floor(GrapeFruit$Exposed*10))*10))

}

\keyword{datasets}
