\name{rmix.norm}
\alias{rmix.norm}
\title{
generate a random sample from a mixture of normals
}
\description{
Generate a random sample from a mixture of normals.
}
\usage{
rmix.norm(n, alpha, mu, sigma = rep(1,length(alpha)))
}
\arguments{
  \item{n}{
sample size.
}
  \item{alpha}{
vector of mixture proportions of length m, the order of the mixture.
}
  \item{mu}{
vector of means of component distributions.
}
  \item{sigma}{
vector of standard deviations of component distributions,
default value: sigma = rep(1,length(alpha)).
}
}
\value{
Return a samples of size n from an m-component normal mixture.
}
\author{
Shaoting Li, Jiahua Chen and Pengfei Li
}
\seealso{
emtest.norm, emtest.norm0, plotmix.norm, plotmix.norm0, pmle.norm, pmle.norm0
}
\examples{
#generate a random sample from a 2 component normal mixture,
#plot the histogram of the data.
x <- rmix.norm(200,c(0.3,0.7),c(-1,2),c(1,2))
hist(x)
}
\keyword{datagen}
\keyword{distribution}