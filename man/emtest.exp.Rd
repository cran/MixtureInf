\name{emtest.exp}
\alias{emtest.exp}
\title{
  Likelihood-based testing for the order of a finite mixture of
  exponentials
}
\description{
  Tests the order of a finite mixture of exponentials with univariate
  mixing parameter.
}
\usage{
emtest.exp(x, m0, C = NULL, len = 5, eps = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    observations whose order is being tested: vector or a matrix with
    the 1st column being the observed values and the 2nd column being
    the corresponding frequency.
  }
  \item{m0}{
    order under null hypothesis
  }
  \item{C}{
    optional tuning parameter for EM-test procedure (if not provided,
    it will be determined by the formulas described in Chen and Li,
    2011)
  }
  \item{len}{
    the number of random restarts of the EM-algorithm
  }
  \item{eps}{
    tolerance value for the convergence of the  EM-algorithm
  }
}
\value{
  Returns an object of class \code{emtest} with the following elements:
  \item{family}{'exponential'}
  \item{m0}{order of the mixture}
  \item{alpha}{estimated mixing proportions}
  \item{theta}{estimated component parameters}
  \item{C}{tuning parameter used}
  \item{ah}{estimated weights in the \eqn{\chi^2} mixture-limiting
    distribution}
  \item{emstat}{EM-test statistics}
  \item{pvalue}{p-values of the EM test statistics}
}

\references{
  Chen, J. and Li, P. (2011). ``Tuning the EM-test for the order of
  finite mixture models.'' \emph{The Canadian Journal of Statistics.}
  39, 389--404.

  Li, P. and Chen, J. (2010). ``Testing the order of a finite mixture
  model.'' \emph{JASA.} 105, 1084--1092.

  Li, P., Chen, J. and Marriott, P. (2009). ``Non-finite Fisher
  information and homogeneity: The EM approach.'' \emph{Biometrika.} 96,
  411--426.
}
\author{
  Jiahua Chen and Pengfei Li.
}
\seealso{
  \code{\link{emtest.binom}}, \code{\link{emtest.norm}},
  \code{\link{emtest.pois}}
}

\examples{
  x <- rexpmix(200, c(0.5, 0.5), c(1, 20))
  emtest.exp(x, 2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
