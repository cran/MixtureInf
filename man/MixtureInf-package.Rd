\name{MixtureInf-package}
\alias{MixtureInf-package}
\alias{MixtureInf}
\docType{package}
\title{
  Likelihood-based testing of the order of a finite mixture
}
\description{
  The \code{MixtureInf} package is an implementation of algorithms for
  likelihood-based testing of the order of a finite mixture, described
  in \cite{Li and Chen (2010)}. For more details refer to the
  bibliography and help for individual functions.
}
\details{
  More information is available on the following topics:
  \itemize{
    \item{\code{emtest.binom}}{test the order of a mixture of binomials}
    \item{\code{emtest.exp}}{test the order of a mixture of exponentials}
    \item{\code{emtest.norm}}{test the order of a mixture of normals
      with variance assumed to be 1}
    \item{\code{emtest.norm2}}{test the order of a mixture of normals
      without equal variance assumption}
    \item{\code{emtest.pois}}{test the order of a mixture of Poissons}
  }
  For a complete list of functions with individual help pages, use:
  \code{library(help=MixtureInf)}.
}
\author{
  Jiahua Chen and Pengfei Li.

  Maintainer: Davor Cubranic \email{cubranic@stat.ubc.ca}
}
\references{
  Chen, J. and Li, P. (2011). ``Tuning the EM-test for the order of
  finite mixture models.'' \emph{The Canadian Journal of Statistics.}
  39, 389--404.

  Li, P. and Chen, J. (2010). ``Testing the order of a finite mixture
  model.'' \emph{JASA.} 105, 1084--1092.

  Li, P. and Chen, J. (2009). ``Non-finite Fisher information and
  homogeneity: The EM approach.'' \emph{Biometrika.} 96, 411--426.

  Chen, J. and Li, P. (2009). ``Hypothesis test for normal mixture
  models: The EM approach.'' \emph{The Annals of Statistics.}
  37, 2523--2542.
}
\keyword{ package }
\keyword{ htest }
\seealso{
  \code{\link{emtest.binom}}, \code{\link{emtest.exp}},
  \code{\link{emtest.norm}}, \code{\link{emtest.norm2}},
  \code{\link{emtest.pois}}. 
}
\examples{
## see examples for function dealing with the distribution of interest
}
