\name{emtest.norm2}
\alias{emtest.norm2}
\title{
  Likelihood-based testing for the order of a finite mixture of normals
}
\description{
  Tests the order of a finite mixture of normals without equal variance assumption
}
\usage{
emtest.norm2(x, m, eps = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    observations whose order is being tested: vector or a matrix with
    the 1st column being the observed values and the 2nd column being
    the corresponding frequency.
  }
  \item{m}{
    order under null hypothesis
  }
  \item{eps}{
    tolerance value for the convergence of the  EM-algorithm
  }
}
\value{
  Returns an object of class \code{emtest} with the following elements:
  \item{family}{'normal_unequalvar'}
  \item{m0}{order of the mixture}
  \item{alpha}{estimated mixing proportions}
  \item{theta}{estimated component parameters \code{mu} and
    \code{sigma}}
  \item{emstat}{EM-test statistics}
  \item{pvalue}{p-values of the EM test statistics}
}

\references{
  Chen, J., Li, P. and Fu, Y. (2012). ``Inference on the order of a
  normal mixture.'' \emph{JASA.}
  107, 1096--1105.

  Chen, J. and Li, P. (2009). ``Hypothesis test for normal mixture
  models: The EM approach.'' \emph{The Annals of Statistics.}
  37, 2523--2542.
}
\author{
  Jiahua Chen and Pengfei Li.
}
\seealso{
  \code{\link{emtest.binom}}, \code{\link{emtest.exp}},
  \code{\link{emtest.norm}}, \code{\link{emtest.pois}}
}

\examples{
  x <- rnormmix(50, c(0.5, 0.5), c(-2, 2), c(1, 4))
  emtest.norm2(x, 2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
