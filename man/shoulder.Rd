\name{Shoulder}

\docType{data}

\alias{Shoulder}

\title{Shoulder flexibility in swimmers}

\description{
This dataset gives the flexibility for the right and left shoulders in 15 swimmers and
15 sedentary people.
}

\usage{data(Shoulder)}

\format{A dataframe with 30 rows and 4 columns.
\tabular{rlll}{
[,1] \tab Subject \tab factor \tab anonymous subjects \cr
[,2] \tab Group \tab factor \tab swimmer or control\cr
[,3] \tab Right \tab numeric \tab right shoulder flexibility (deg.)\cr
[,4] \tab Left \tab numeric \tab left shoulder flexibility (deg.)\cr
}
}

\source{Private communication. Karine Monteil, CRIS, Lyon 1 University, FRANCE.}

\references{Monteil, K., Taiar, R., Champely, S. and Martin, J. (2002) Competitive swimmers versus sedentary people: 
a predictive model based upon normal shoulders flexibility. Journal of Human Movement Studies, 43 , 17-34.
}

\examples{
data(Shoulder)

# Is there some heteroscedasticity?
with(Shoulder,plot(paired(Left,Right)))

# Swimmers are indeed quite different
with(Shoulder,plot(paired(Right,Left),groups=Group))

# A first derived variable to compare the amplitude in flexibilty
with(Shoulder,boxplot(((Left+Right)/2)~Group,ylab="mean shoulder flexibility"))

# A second derived variable to study shoulder asymmetry
with(Shoulder,boxplot((abs(Left-Right))~Group,ylab="asymmetry in shoulder flexibility"))
}

\keyword{datasets}
