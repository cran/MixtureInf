\name{pearson}
\docType{data}
\alias{pearson}
\title{
Pearson's crab data
}
\description{
This data set contains the ratio of "forehead" breadth to body length for 1000 crabs sampled 
at Naples by Professor W.F.R. Weldon. 
Pearson (1894) used a two component normal mixture model to fit this data set.
The data frame has 29 rows and 2 columns.
}
\usage{
#Load this data set for subsequent analysis:
data(pearson)
}
\format{
This data frame contains the following columns:

ratio: the boundaries of grouping intervals.

freq: the frequencies of observation falling into each interval.
}
\references{
Pearson K (1894). Contributions to the mathematical theory of evolution. Philosophical Transactions of the Royal Society of London. A, 185, 71--110.
}
\examples{
#fit a 2 component normal mixture model to this data set.
data(pearson)
pmle.norm(pearson,2,1)
}
\keyword{datasets}