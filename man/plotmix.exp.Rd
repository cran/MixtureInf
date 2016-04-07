\name{plotmix.exp}
\alias{plotmix.exp}
\title{
Histogram of the observations and the fitted density of a mixture of exponentials
}
\description{
Plot the histogram of the observations and the fitted density of a mixture of exponentials.
}
\usage{
plotmix.exp(x, theta, hist = 1, comp = TRUE, k = 20, h = 1, main = "",
            xlab = "Observations", ylab = "")
}
\arguments{
  \item{x}{
data, can be either a vector or a matrix with the 1st column being the observed values and the 2nd column being the corresponding frequencies.
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
Return the histogram of observations and the plot of the fitted density.
}
\author{
Shaoting Li, Jiahua Chen and Pengfei Li
}
\seealso{
emtest.exp, pmle,exp, rmix.exp
}
\examples{
#load the lengthofstay data set,
#fit a 2 compoent exponential mixture model, 
#plot the histogram of the observations and the fitted density.
data(lengthofstay)
out <- pmle.exp(lengthofstay,2)
plotmix.exp(lengthofstay,out)
}
\keyword{dplot}
