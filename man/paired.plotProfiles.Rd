\name{paired.plotProfiles}
\alias{paired.plotProfiles}

\title{Profile plot}
\description{
Produce a profile plot or before-after plot  or 1-1 plot for paired data.
}
\usage{
paired.plotProfiles(df, condition1, condition2, groups = NULL,subjects, 
    facet = TRUE, ...) 
}

\arguments{
  \item{df}{
a data frame.
}
  \item{condition1}{
name of the variable corresponding to the second sample.
}
  \item{condition2}{
name of the variable corresponding to the first sample.
}

  \item{groups}{
names of the variable corresponding to groups (optional).
}
  \item{subjects}{
names of the variable corresponding to subjects.
}
  \item{facet}{
faceting or grouping strategy for plotting?
}

\item{\dots}{
further arguments to be passed to methods.}
}


\value{
a graphical object of class ggplot.
}


\references{Cox, N.J. (2004) Speaking data: graphing agreement and disagreement. The Stata Journal, 4, 329-349.}


\author{Stephane CHAMPELY}

\seealso{plotBA,plotMcNeil}

\examples{
data(PrisonStress)
paired.plotProfiles(PrisonStress,"PSSbefore","PSSafter",subjects="Subject",groups="Group")

# Changing the line colour
paired.plotProfiles(PrisonStress,"PSSbefore","PSSafter")+geom_line(colour="red")
}