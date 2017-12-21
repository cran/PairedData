\name{GDO}

\docType{data}

\alias{GDO}

\title{Agreement study}

\description{
This dataset gives the same measurements of muscle activation (EMG) in 3 days corresponding to a reproductibility study
for 18 tennis players.
}

\usage{data(GDO)}

\format{A dataframe with 18 rows and 4 columns.
\tabular{rlll}{
[,1] \tab Subject \tab factor \tab anonymous subjects \cr
[,2] \tab Day1 \tab numeric \tab measurement first day\cr
[,3] \tab Day2 \tab numeric \tab measurement second day\cr
[,4] \tab Day3 \tab numeric \tab measurement third day\cr
}
}

\source{Private communication. Samuel Rota, CRIS, Lyon 1 University, FRANCE}

\seealso{
packages: agreement, irr and MethComp.
}
\examples{
data(GDO)

# Building new vectors for performing
# a repeated measures ANOVA
# with a fixed Day effect
Activation<-c(GDO[,2],GDO[,3],GDO[,4])
Subject<-factor(rep(GDO[,1],3))
Day<-factor(rep(c("D1","D2","D3"),rep(18,3)))
aovGDO<-aov(Activation~Day+Error(Subject))
summary(aovGDO)

# Reliability measurement: SEM and ICC(3,1)
sqrt(12426)
72704/(72704+12426)
}

\keyword{datasets}
