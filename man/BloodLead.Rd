\name{BloodLead}

\docType{data}

\alias{BloodLead}

\title{Blood lead levels data from Pruzek & Helmreich (2009)}

\description{
This dataset presents matched paired data corresponding to blood lead levels for 33 children of parents who had worked in a lead related factory and 33 control children from their neighborhood.
The two samples have different dispersions and their correlation is small. }

\usage{data(BloodLead)}

\format{A dataframe with 33 rows and 3 columns: 
\tabular{rlll}{
[,1] \tab Pair \tab factor \tab matched pair of chidren \cr
[,2] \tab Exposed \tab numeric \tab blood lead levels (mg/dl) for exposed children \cr
[,3] \tab Control \tab numeric \tab blood lead levels for controls \cr
}
}

\source{
Morton, D., Saah, A., Silberg, S., Owens, W., Roberts, M. & Saah, M. (1982) Lead absorption in children of employees in a lead related industry. American Journal of Epimediology, 115, 549-55.
}

\references{
Pruzek, R.M. & Helmreich, J.E. (2009) Enhancing dependent sample analysis with graphics. Journal of Statistics Education, 17 (1).
}

\examples{
data(BloodLead)

# Control values are clealy less dispersed (and inferior)
# than exposed levels
with(BloodLead,plot(paired(Control,Exposed),type="McNeil"))
with(BloodLead,Var.test(paired(Control,Exposed)))

with(BloodLead,grambsch.Var.test(paired(Control,Exposed)))
with(BloodLead,bonettseier.Var.test(paired(Control,Exposed)))

# Correlation is small (bad matching)
with(BloodLead,cor.test(Control,Exposed))
with(BloodLead,winsor.cor.test(Control,Exposed))
}


\keyword{datasets}
