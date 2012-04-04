\name{paired-class}
\Rdversion{1.1}
\docType{class}
\alias{paired-class}


\title{Class \code{"paired"}}
\description{
An object of class paired is a dataframe with two columns sharing the same class (usually numeric).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("paired", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
    \item{\code{names}:}{Object of class \code{"character"} ~~ }
    \item{\code{row.names}:}{Object of class \code{"data.frameRowLabels"} ~~ }
    \item{\code{.S3Class}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{data.frame}"}, directly.
Class \code{"\linkS4class{list}"}, by class "data.frame", distance 2.
Class \code{"\linkS4class{oldClass}"}, by class "data.frame", distance 2.
Class \code{"\linkS4class{vector}"}, by class "data.frame", distance 3.
}
\section{Methods}{
  \describe{
    \item{effect.size}{\code{signature(object = "paired")}: ... }
    \item{summary}{\code{signature(object = "paired")}: ... }
    \item{plot}{\code{signature(object = "paired")}: ... }

	 }
}

\author{
Stephane Champely
}
\examples{
data(iceskating)
p<-with(iceskating,paired(Extension,Flexion))
show(p)
plot(p)
summary(p)
}
\keyword{classes}
