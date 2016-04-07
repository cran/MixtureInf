\name{rmix.exp}
\alias{rmix.exp}
\title{
generate a random sample from a mixture of exponentials
}
\description{
Generate a random sample from a mixture of exponentials.
}
\usage{
rmix.exp(n, alpha, theta)
}
\arguments{
  \item{n}{
sample size
}
  \item{alpha}{
vector of mixture proportions, with length of m, the order of the mixture.
}
  \item{theta}{
vector of parameters of component distributions.
}
}
\value{
Return a sample of size n from an m-component exponential mixture.
}
\author{
Shaoting Li, Jiahua Chen and Pengfei Li
}
\seealso{
emtest.exp, plotmix.exp, pmle.exp
}
\examples{
#generate a random sample from a 2 component exponential mixture,
#compute the sample mean and variance.
x <- rmix.exp(200,c(0.3,0.7),c(2,8))
mean(x)
var(x)
}
\keyword{datagen}
\keyword{distribution}
