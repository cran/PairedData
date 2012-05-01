\name{paired.plotCor}
\alias{paired.plotCor}

\title{Paired correlation plot}
\description{
Produce a squared scatterplot for paired data (same units for both axes), including the first bisector line for reference.
}
\usage{
paired.plotCor(df, condition1, condition2, groups = NULL, 
facet = TRUE, ...)
}

\arguments{
  \item{df}{
a data.frame.
}
  \item{condition1}{
name of the variable corresponding to the first sample.
}
  \item{condition2}{
name of the variable corresponding to the first sample.
}
\item{groups}{
name of the variable corresponding to the groups (optional).
}

  \item{facet}{
faceting or grouping strategy for plotting?
}

\item{\dots}{
arguments to be passed to methods}
}

\value{
a graphical object of class ggplot.
}


\author{Stephane CHAMPELY}


\examples{
data(PrisonStress)
paired.plotCor(PrisonStress,"PSSbefore","PSSafter")

# Changing the theme of the ggplot object
paired.plotCor(PrisonStress,"PSSbefore","PSSafter")+theme_bw()
}

