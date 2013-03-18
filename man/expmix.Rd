\name{rexpmix}
\alias{rexpmix}

\title{Random sample from a mixture of exponentials}
\description{
  Generate a random sample for a mixture of univariate exponential
  distributions.
}
\usage{
rexpmix(n, alpha, scale)
}

\arguments{
  \item{n}{number of samples}
  \item{alpha}{vector of mixture probabilities, with length equal to
    \eqn{m}, the desired number of components.}
  \item{scale}{vector of scales (inverse of rate) of each component.}
}

\value{
  Returns a vector of length \eqn{n} of samples from an
  \eqn{m}-component mixture of univariate exponential distributions.
}

\seealso{
  Functions for mixtures of other standard distributions, such as
  \code{\link{rnormmix}} for the normal and \code{\link{rpoismix}}
  for the Poisson distribution.
}

\examples{
# Draw random samples from a mixture of two equally likely exponential
# distributions, one with scale = 1, the other with scale = 0.05
x <- rexpmix(200, c(0.5, 0.5), c(1, 20))
hist(x, probability=TRUE)
lines(density(x), col='red')
}

\keyword{ distribution }
