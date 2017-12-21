\name{Tobacco}

\docType{data}

\alias{Tobacco}

\title{Tobacco data from Snedecor and Cochran (1967)}

\description{
This dataset presents 8 paired data corresponding to numbers of lesions caused by two virus preparations inoculated into the two halves of each tobacco leaves.
}

\usage{data(Tobacco)}

\format{A dataframe with 8 rows and 3 columns:  
\tabular{rlll}{
[,1] \tab Plant \tab factor \tab \cr
[,2] \tab Preparation_1 \tab numeric \tab number of lesions \cr
[,3] \tab Preparation_2 \tab numeric \tab number of lesions \cr
}
}

\source{
Snedecor, G.W. and Cochran, W.G. (1967) Statistical Methods, 6th ed. Iowa State University Press: Ames.
}

\references{
\itemize{
\item Pruzek, R.M. & Helmreich, J.E. (2009) Enhancing dependent sample analysis with graphics. Journal of Statistics Education, 17 (1).
\item Preece, D.A. (1982) t is for trouble (and textbooks): a critique of some examples of the paired-samples t-test. The Statistician, 31 (2), 169-195.
}
}

\examples{
data(Tobacco)

# A clear outlier
with(Tobacco,plot(paired(Preparation_1,Preparation_2)))


# Comparison of normal and robust tests
with(Tobacco,t.test(paired(Preparation_1,Preparation_2)))
with(Tobacco,yuen.t.test(paired(Preparation_1,Preparation_2)))

with(Tobacco,Var.test(paired(Preparation_1,Preparation_2)))
with(Tobacco,grambsch.Var.test(paired(Preparation_1,Preparation_2)))

with(Tobacco,cor.test(Preparation_1,Preparation_2))
with(Tobacco,winsor.cor.test(Preparation_1,Preparation_2))

# Maybe a transformation
require(MASS)
with(Tobacco,eqscplot(log(Preparation_1),log(Preparation_2)))
abline(0,1,col="red")
}

\keyword{datasets}
