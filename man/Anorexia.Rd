\name{Anorexia}

\docType{data}

\alias{Anorexia}

\title{Anorexia data from Pruzek & Helmreich (2009)}

\description{
This dataset presents 17 paired data corresponding to the weights of girls before and after treatment for anorexia. A more complete version can be found in the package MASS.
There is actually a cluster of four points in this dataset.
}

\usage{data(Anorexia)}

\format{A dataframe with 17 rows and 2 numeric columns: 
\tabular{rlll}{
[,1] \tab Prior \tab numeric \tab weight (lbs) before therapy \cr
[,2] \tab Post \tab numeric \tab weight (lbs) after therapy \cr
}
}

\seealso{anorexia in MASS}

\source{
Hand, D.J., McConway, K., Lunn, D. & Ostrowki, editors (1993) A Handbook of Small Data Sets. Number 232, 285. Chapman & Hall: New-York.
}

\references{
Pruzek & Helmreich (2009) Enhancing dependent sample analysis with graphics. Journal of Statistics Education, 17 (1).
}

\examples{
data(Anorexia)

# Visualization of the cluster
with(Anorexia,plot(paired(Prior,Post),type="profile"))

# The effects of trimming or winsorizing 
# with 4 outliers (n=17)
17*0.2
with(Anorexia,summary(paired(Prior,Post)))
17*0.25
with(Anorexia,summary(paired(Prior,Post),tr=0.25))
}


\keyword{datasets}
