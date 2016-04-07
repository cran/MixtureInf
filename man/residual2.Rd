\name{residual2}
\docType{data}
\alias{residual2}
\title{
residual diagnostics for mixture models
}
\description{
This data set is from Lindsay and Roeder (1992). 
It includes the number of boys in families of size 12.
The number of families is 6115. The data frame has 13 rows and 2 columns.
}
\usage{
#Load this data set for subsequent analysis:
data(residual2)
}
\format{
This data frame contains the 2 columns, 

count: number of boys in family.

freq: number of families with corresponding.
}
\references{
Lindsay, B. G. and Roeder, K. (1992). Residual diagnostics for mixture models.
Journal of the American Statistical Association,87(419), 785--794.
}
\examples{
#fit a 2 component binomial mixture model to this data set.
data(residual2)
pmle.binom(residual2,12,2,1)
}
\keyword{datasets}