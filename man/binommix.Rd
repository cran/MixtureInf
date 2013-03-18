\name{rbinommix}
\alias{rbinommix}

\title{Random sample from a mixture of binomials}
\description{
  Generate a random sample for a mixture of univariate binomial
  distributions.
}
\usage{
rbinommix(n, alpha, prob, size)
}

\arguments{
  \item{n}{number of samples}
  \item{alpha}{vector of mixture probabilities, with length equal to
    \eqn{m}, the desired number of components.}
  \item{prob}{vector of probabilities of success of each component.}
  \item{size}{number of trials}
}

\value{
  Returns a vector of length \eqn{n} of samples from an
  \eqn{m}-component mixture of univariate binomial distributions.
}

\seealso{
  Functions for mixtures of other standard distributions, such as
  \code{\link{rnormmix}} for the normal and \code{\link{rpoismix}}
  for the Poisson distribution.
}

\examples{
# Draw random samples from a mixture of two equally likely binomials,
# one with p(success) = 0.1, the other with p(success) = 0.9
x <- rbinommix(200, c(0.5, 0.5), c(0.1, 0.9), 16)
hist(x, probability=TRUE)
lines(density(x), col='red')
}

\keyword{ distribution }
