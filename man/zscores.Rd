\name{zscores}
\docType{data}
\alias{zscores}
\title{
z-scores
}
\description{
This data set consists of 6033 z-scores, transformed from the two-sample t-test statistics
based on the prostate cancer data set of Singh et al. (2002).
The prostate cancer data set consists of gene expression levels of 
6033 genes of 52 prostate cancer patients and 50 normal control subjects.
By Efron (2004), if the gene is not differentially expressed, the corresponding z should follow N(0, 1). 
Normal mixture with known variance 1 is suitable to this data set. 
The data frame has 6033 rows and 1 columns.
}
\usage{
#Load this data set for subsequent analysis:
data(zscores)
}
\format{
This data frame contains one column:

z-scores: 6033 z-scores.
}
\references{
Efron, B. (2004), Large-Scale Simultaneous Hypothesis Testing: The Choice of a
Null Hypothesis, Journal of the American Statistical Association, 99, 96--104.

Singh, D., Febbo, P. G., Ross, K., Jackson, D. G., Manola, J., Ladd, C., Tamayo,
P., Renshaw, A. A., D'Amico, A. V., Richie, J. P., Lander, E. S., Loda,
M., Kantoff, P. W., Golub, T. R., and Sellers, W. R. (2002), Gene Expression
Correlates of Clinical Prostate Cancer Behavior, Cancer Cell, 1, 203--209.
}
\examples{
\dontrun{
#fit a mixture of 2 normals with equal and known variance 1 to this data set.
data(zscores) 
pmle.norm0(zscores,var=1,2,1)
}
}
\keyword{datasets}
