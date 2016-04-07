\name{agesch}
\docType{data}
\alias{agesch}
\title{
age of onset of schizophrenia
}
\description{
This data set is from a schizophrenia study by Levine (1981), who
collacted the results of seven studies on the age of onset of schizophrenia of 99
females and 152 males. This data set contains observed values of only male data.
Everitt et al., (2001) fitted the male data using a two-component normal mixture model. 
The data frame has 152 rows and 1 column.
}
\usage{
#Load this data set for subsequent analysis:
data(agesch)
}
\format{
This data frame contains one column:

age: age of onset of schizophrenia of 152 males.
}
\references{
Levine, R. (1981). Sex differences in schizophrenia: Timing or subtypes? Psychological Bulletin, 90, 432--444.

Everitt, B. S., Landau, S., and Leese, M. (2001). Cluster Analysis (4th ed.). New York, NY: Oxford University Press, Inc.
}
\examples{
#fit a finite normal mixture model of order 2 to this data set,
#plot the histogram of the observations and the fitted density function.
data(agesch)
out <- pmle.norm(agesch,2,1)
plotmix.norm(agesch,out)
}
\keyword{datasets}