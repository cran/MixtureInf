\name{timesoffailure}
\docType{data}
\alias{timesoffailure}
\title{
successive failure times
}
\description{
The data set contains the times of successive failures for the air conditioning system
of each member in a fleet of 13 Boeing 720 jet aircrafts. The data frame has 213 rows
and 1 column.
}
\usage{
#Load this data set for subsequent analysis:
data(timesoffailure)
}
\format{
This data frame contains one column,

x: pooled failure times about 213 observations.
}
\references{
Proschan, F. (1963). Theoretical explanation of observed decreasing failure rate. Technometrics 5, 375--83.
}
\examples{
#fit a finite exponential mixture model of order 2 to this data set. 
data(timesoffailure)
pmle.exp(timesoffailure,2,1)
}
\keyword{datasets}