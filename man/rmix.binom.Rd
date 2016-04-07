\name{rmix.binom}
\alias{rmix.binom}
\title{
generate a random sample from a mixture of binomials
}
\description{
Generate a random sample from a mixture of binomials.
}
\usage{
rmix.binom(n, alpha, theta, size)
}
\arguments{
  \item{n}{
sample size
}
  \item{alpha}{
vector of mixture proportions of length m, the order of the mixture.
}
  \item{theta}{
vector of probabilities of success of component distributions.
}
  \item{size}{
number of trials.
}
}
\value{
Return sample of size n from an m-component binomial mixture.
}
\author{
Shaoting Li, Jiahua Chen and Pengfei Li
}
\seealso{
emtest.binom, plotmix.binom, pmle.binom
}
\examples{
#generate a random sample from a 2 component binomial mixture,
#plot the histogram of the data
x <- rmix.binom(200,c(0.3,0.7),c(0.4,0.6),20)
hist(x)
}
\keyword{datagen}
\keyword{distribution}