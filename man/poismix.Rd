\name{rpoismix}
\alias{rpoismix}

\title{Random sample from a mixture of Poissons}
\description{
  Generate a random sample for a mixture of univariate Poisson
  distributions.
}
\usage{
rpoismix(n, alpha, lambda)
}

\arguments{
  \item{n}{number of samples}
  \item{alpha}{vector of mixture probabilities, with length equal to
    \eqn{m}, the desired number of components.}
  \item{lambda}{vector of means of each component.}
}

\value{
  Returns a vector of length \eqn{n} of samples from an
  \eqn{m}-component mixture of univariate Poisson distributions.
}

\seealso{
  Functions for mixtures of other standard distributions, such as
  \code{\link{rnormmix}} for the normal and \code{\link{rexpmix}}
  for the exponential distribution.
}

\examples{
# Draw random samples from a mixture of two equally likely Poisson
# distributions, one with mean = 1, and the other with mean = 10
x <- rpoismix(200, c(0.5, 0.5), c(1,10))
hist(x, probability=TRUE)
lines(density(x), col='red')
}

\keyword{ distribution }
