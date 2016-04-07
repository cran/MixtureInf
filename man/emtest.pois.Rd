\name{emtest.pois}
\alias{emtest.pois}
\title{
test the order of a mixture of Poissons
}
\description{
EM-test for the order of a finite mixture of Poissons.
}
\usage{
emtest.pois(x, m0 = 1, C = NULL, 
	    inival = NULL, len = 10, niter = 50, tol = 1e-6, k = 3, rformat = FALSE)
}
\arguments{
  \item{x}{
data, can be either a vector or a matrix with the 1st column being the observed values and the 2nd column being the corresponding frequencies.
}
  \item{m0}{
order of the finite mixture model under the null hypothesis, default value: m0 = 1.
}
  \item{C}{
optional tuning parameter for EM-test procedure, default value: C = NULL.
(if not provided, it will be determined by the formulas described in Chen and Li, 2011).
}
  \item{inival}{
initial values for the EM-algorithm to compute the MLE under the null model, 
a 2m0-dimension vector including m0 mixing proportions and m0 component parameters, 
or a matrix with 2m0 columns, default value: inival = NULL. 
(if not provided, random initial values are used.)
}
  \item{len}{
number of random initial values for the EM-algorithm, default value: len = 10.
}
  \item{niter}{
number of iterations for all initial values in the EM-algorithm. 
The algorithm runs EM-iteration niter times from each initial value. 
The iteration will restart from the parameter value with the highest 
likelihood value at the point and run until convergence. default value: niter = 50.
}
  \item{tol}{
tolerance level for the convergence of the EM-algorithm, default value: tol = 1e-6.
}
  \item{k}{
number of EM iterations: default value: k = 3.
}
  \item{rformat}{
form of the digital output: default of R package is used when rformat = T; If rformat = T, the digital 
output is rounded to the 3rd dicimal place if it is larger than 0.001, keeps 3 significant digits otherwise. 
The default value of rformat is F.
}
}
\value{
Return an object of class EM-test with the following elements:

MLE of the parameters under the null hypothesis (order = m0)

Parameter estimates under the specific alternative whose order is 2m0.

EM-test statistic

P-value

Level of penalty
}
\references{
Chen, J. and Li, P. (2011). Tuning the EM-test for the order of finite mixture models. The Canadian Journal of Statistics. 39, 389--404.

Li, P. and Chen, J. (2010). Testing the order of a finite mixture model. JASA. 105, 1084--1092.

Li, P., Chen, J. and Marriott, P. (2009). Non-finite Fisher information and homogeneity: The EM approach. Biometrika. 96, 411--426.
}
\author{
Shaoting Li, Jiahua Chen and Pengfei Li
}
\seealso{
plotmix.pois, pmle.pois, rmix.pois
}
\examples{
#load the earthquake data set,
#conduct homogeneity test by the EM-test.
data(earthquake)
emtest.pois(earthquake,2)
}
\keyword{htest}
