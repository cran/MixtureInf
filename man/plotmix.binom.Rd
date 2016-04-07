\name{plotmix.binom}
\alias{plotmix.binom}
\title{
Histogram of the observations and the fitted probability mass function of a mixture of binomials
}
\description{
Plot the histogram of the observations and the fitted probability mass function of a mixture of binomials.
}
\usage{
plotmix.binom(x, size, theta, hist = 1, comp = TRUE, h = 1, main = "",
              xlab = "Observations", ylab = "")
}
\arguments{
  \item{x}{
data, can be either a vector or a matrix with the 1st column being the observed values and the 2nd column being the corresponding frequencies.
}
  \item{size}{
number of trials.
}
  \item{theta}{
parameter values, output of emtest.binom or pmle.binom, or a vector of parameter values, include mixing proportions and component parameters.
}
  \item{hist}{
style of histogram, hist = 0,1,2 are permitted. hist = 0 means no histogram, hist = 1,2 give two styles, default value: hist = 1.
}
  \item{comp}{
a parameter specifies the form of the component probability functions in the graph.
The component probability functions will be plotted only if the order is 3 or lower.
comp = T switches on the component probability functions in the graph. If comp=F, no component probability functions will be included. 
default value: comp = T.
}
  \item{h}{
a number between 0 and 1, specifying the height of the probability function, default value: h = 1 (complete probability function).
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
Return the histogram of the observations and the plot of the fitted probability mass function.
}
\author{
Shaoting Li, Jiahua Chen and Pengfei Li
}
\seealso{
emtest.binom, pmle.binom, rmix.binom
}
\examples{
#generate a random sample from a 2 component binomial mixture,
#compute the PMLE of parameters based on the 2 component binomial mixture model,
#plot the histgoram of the observations and the fitted probability mass function.
x <- rmix.binom(200,c(0.3,0.7),c(0.4,0.6),50)
out <- pmle.binom(x,50,2)
plotmix.binom(x,50,out)
}
\keyword{dplot}
