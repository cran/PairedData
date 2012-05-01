\name{paired}
\alias{paired}

\title{
Paired
}
\description{
This function creates objects of class paired
}
\usage{
paired(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
first vector.
}
  \item{y}{
second vector.
}
}
\details{
The two vectors must share the same class. Moreover, for vectors of class factor, they must have the same levels.
}
\value{
An object of class paired.
}
\author{
Stephane Champely
}
\examples{
x<-rnorm(15)
y<-rnorm(15)
p1<-paired(x,y)
show(p1)

data(IceSkating)
p2<-with(IceSkating,paired(Extension,Flexion))
show(p2)

}
