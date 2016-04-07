\name{plotmix.norm0}
\alias{plotmix.norm0}
\title{
Histogram of the observations and the fitted density of a mixture of normals with equal and known variance
}
\description{
Plot the histogram of the observations and the fitted density of a mixture of normals with equal and known variance.
}
\usage{
plotmix.norm0(x, var, theta, hist = 1, comp = TRUE, k = 20, h = 1, main = "",
              xlab = "Observations", ylab = "")
}
\arguments{
  \item{x}{
data, can be either a vector or a matrix with the 1st column being the observed values and the 2nd column being the corresponding frequencies.
}
  \item{var}{
known component variance.
}
  \item{theta}{
parameter values, output of emtest.exp or pmle.exp, or a vector of parameter values, include mixing proportions and component parameters.
}
  \item{hist}{
style of histogram. hist = 0,1,2 are permitted. hist = 0 means no histogram, hist = 1,2 give two styles, default value: hist = 1.
}
  \item{comp}{
a parameter specifies the form of the component densities in the graph.
comp = T switches on the component densities in the graph. If comp = F, no component densities will be included.
default value: comp = T.
}
  \item{k}{
number of bins for the histogram, default value: k = 20.
}
  \item{h}{
a number between 0 and 1, specifying the height of the density function, default value: h = 1 (complete density function).
}
  \item{main}{
title of graph, default value: main = "".
}
  \item{xlab}{
label of x-axis, default value: xlab = "Observations".
}
  \item{ylab}{
label of y-axis, default value: ylab = "".
}
}
\value{
Return the histogram of the observations and the plot of the fitted density.
}
\author{
Shaoting Li, Jiahua Chen and Pengfei Li
}
\seealso{
emtest.norm0, pmle.norm0, rmix.norm
}
\examples{
#generate a random sample from a 2 component normal mixture,
#compute the PMLE of parameters under a 2 component normal mixture model,
#plot the histgoram of the observations and the fitted density function.
x <- rmix.norm(200,c(0.3,0.7),c(-2,3))
out <- pmle.norm0(x,var=1,2)
plotmix.norm0(x,var=1,out)
}
\keyword{dplot}
