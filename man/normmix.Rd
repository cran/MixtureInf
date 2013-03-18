\name{rnormmix}
\alias{rnormmix}

\title{Random sample from a mixture of normals}
\description{
  Generate a random sample for a mixture of univariate normal
  distributions.
}
\usage{
rnormmix(n, alpha, mean, sd = rep(1, length(alpha)))
}

\arguments{
  \item{n}{number of samples}
  \item{alpha}{vector of mixture probabilities, with length equal to
    \eqn{m}, the desired number of components.}
  \item{mean}{vector of means of each component.}
  \item{sd}{vector of standard deviations of each component.}
}

\value{
  Returns a vector of length \eqn{n} of samples from an
  \eqn{m}-component mixture of univariate normal distributions.
}

\seealso{
  Functions for mixtures of other standard distributions, such as
  \code{\link{rbinommix}} for the binomial and \code{\link{rexpmix}}
  for the exponential distribution.
}

\examples{
# Draw random samples from a mixture of two equally likely normal
# distributions, one with mean = -2, s.d. = 1, and the other with
# mean=2, s.d. = 4
x <- rnormmix(200, c(0.5, 0.5), c(-2, 2), c(1, 4))
hist(x, probability=TRUE)
lines(density(x), col='red')
}

\keyword{ distribution }
