\name{PrisonStress}

\docType{data}

\alias{PrisonStress}

\title{Stress in prison}

\description{
This dataset gives the PSS (stress measurement) for 26 people in prison at the entry and at the exit. 
Part of these people were physically trained during their imprisonment.
}

\usage{data(PrisonStress)}

\format{A dataframe with 26 rows and 4 columns.
\tabular{rlll}{
[,1] \tab Subject \tab factor \tab anonymous subjects \cr
[,2] \tab Group \tab factor \tab sport or control\cr
[,3] \tab PSSbefore \tab numeric \tab stress measurement before training\cr
[,4] \tab PSSafter \tab numeric \tab stress measurement after training\cr
}
}

\source{Private communication. Charlotte Verdot, CRIS, Lyon 1 University, FRANCE}

\references{
Verdot, C., Champely, S., Massarelli, R. and Clement, M. (2008) Physical activities in prison as 
a tool to ameliorate detainees mood and well-being. 
International Review on Sport and Violence, 2.
}

\examples{
data(PrisonStress)

# The two groups are not randomized! 
# The control group is less stressed before the experiment
with(PrisonStress,boxplot(PSSbefore~Group,ylab="Stress at the eginning of the study"))

# But more stressed at the end!
with(PrisonStress,boxplot(PSSafter~Group,ylab="22 weeks later"))

# So the effects of physical training seems promising
with(PrisonStress,plot(paired(PSSbefore,PSSafter),groups=Group,type="BA",facet=FALSE))

# Testing using gain scores analysis
difference<-PrisonStress$PSSafter-PrisonStress$PSSbefore
t.test(difference~PrisonStress$Group,var.equal=TRUE)

# Testing using ANCOVA
lmJail<-lm(PSSafter~PSSbefore*Group,data=PrisonStress)
anova(lmJail)

# Testing using repeated measures ANOVA
PSS<-c(PrisonStress$PSSbefore,PrisonStress$PSSafter)
Time<-factor(rep(c("Before","After"),c(26,26)))
Subject<-rep(PrisonStress$Subject,2)
Condition<-rep(PrisonStress$Group,2)
aovJail<-aov(PSS~Condition*Time+Error(Subject))
summary(aovJail)
}

\keyword{datasets}
