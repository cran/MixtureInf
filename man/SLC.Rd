\name{SLC}
\docType{data}
\alias{SLC}
\title{
SLC (sodium-lithium countertransport)
}
\description{
This data set contains 190 SLC measurements studied by Roeder, K. (1994). 
Roeder (1994) analyzed this data and concluded that a three component normal mixture
with equal variance is most suitable. 
Chen et al. (2012) also analyzed this data and gave a slightly better fit by 
a two component normal mixture with unequal variance.
The data frame has 190 rows and 1 column.
}
\usage{
#Load this data set for subsequent analysis:
data(SLC)
}
\format{
This data frame contains one column:

SLC: 190 slc measurements.
}
\references{
Chen, J., Li, P. and Fu, Y. (2012). Inference on the order of a normal mixture. JASA. 107, 1096--1105.

Roeder, K. (1994), A Graphical Technique for Determining the Number of Components
in a Mixture of Normals, Journal of the American Statistical Association, 89,
487--500.
}
\examples{
#fit a 2 component normal mixture model to this data set, 
#test for order = 2 using EM-test.
data(SLC)
a <- c(0.6,0.4,0.2,0.3,0.01,0.01)
pmle.norm(SLC,2,1,inival=a)
emtest.norm(SLC,2,inival=a,len=5)
}
\keyword{datasets}
