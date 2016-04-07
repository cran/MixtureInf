\name{rmix.pois}
\alias{rmix.pois}
\title{
generate a random sample from a mixture of Possions
}
\description{
Generate a random sample from a mixture of Possions.
}
\usage{
rmix.pois(n, alpha, theta)
}
\arguments{
  \item{n}{
sample size.
}
  \item{alpha}{
vector of mixture proportions of length m, the order of the mixture.
}
  \item{theta}{
vector of parameters of each component.
}
}
\value{
Return a sample of size n from an m-component Poisson mixture.
}
\author{
Shaoting Li, Jiahua Chen and Pengfei Li
}
\seealso{
rmix.binom, rmix.exp, rmix.norm
}
\examples{
#generate a random sample from a 2 component Poisson mixture,
#and compute the sample mean and variance.
x <- rmix.pois(200,c(0.3,0.7),c(2,5))
mean(x)
var(x)
}
\keyword{datagen}
\keyword{distribution}
