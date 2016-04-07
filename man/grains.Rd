\name{grains}
\docType{data}
\alias{grains}
\title{
grains data
}
\description{
The data set contains the square root of the total number of grains for each planty from Loisel et al., (1994). 
Loisel et al., (1994) suggested using a finite normal mixture model.
The grains data frame has 150 rows and 1 column.
}
\usage{
#Load this data set for subsequent analysis:
data(grains)
}
\format{
This data frame contains one column:

x: square root of the total number of grains for each planty.
}
\references{
Loisel, P., Goffinet, B., Monod, H., and Montes De Oca, G. (1994). Detecting a major gene in an F2 population. Biometrics, 50, 512--516.
}
\examples{
#fit a finite noraml mixture model of order 2 to this data set. 
data(grains)
pmle.norm(grains,2,1)
}
\keyword{datasets}