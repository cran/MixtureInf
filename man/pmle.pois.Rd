\name{pmle.pois}
\alias{pmle.pois}
\title{
compute the PMLE or MLE of the parameters under a mixture of Possions
}
\description{
Compute the PMLE or MLE of the parameters under a mixture of Possions. When the level of penalty is 0, PMLE reduces to MLE.
}
\usage{
pmle.pois(x, m0 = 1, lambda = 0, 
	  inival = NULL, len = 10, niter = 50, tol = 1e-06, rformat = FALSE)
}
\arguments{
  \item{x}{
data, can be either a vector or a matrix with the 1st column being the observed values and the 2nd column being the corresponding frequencies.
}
  \item{m0}{
order of the finite mixture model, default value: m0 = 1.
}
  \item{lambda}{
level of penalty, default value: lambda = 0.
}
  \item{inival}{
initial values for the EM-algorithm, a 2m0-dimension vector 
including m0 mixing proportions and m0 component parameters, 
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
  \item{rformat}{
form of the digital output: default of R package is used when rformat = T; If rformat = T, the digital 
output is rounded to the 3rd dicimal place if it is larger than 0.001, keeps 3 significant digits otherwise. 
The default value of rformat is F.
}
}
\value{
Return the PMLE or MLE of the parameters with order = m0 (mixing proportions and component parameters), 
log-likelihood value at the PMLE or MLE and the penalized log-likelihood value at the PMLE.
}
\author{
Shaoting Li, Jiahua Chen and Pengfei Li
}
\seealso{
emtest.pois, pmle.pois, rmix.pois
}
\examples{
#generate a random sample from a 2 component Poisson mixture,
#compute the PMLE of the parameters under the 2 component Poisson mixture model,
#plot the histgoram of the observations and the fitted probability mass function.
x <- rmix.pois(200,c(0.3,0.7),c(2,5))
out <- pmle.pois(x,2)
plotmix.pois(x,out)
}
\keyword{models}
