\name{plotmix.pois}
\alias{plotmix.pois}
\title{
Histogram of the observations and the fitted probability mass function of a mixture of Poissons
}
\description{
Plot the histogram of the observations and the fitted probability mass function of a mixture of Poissons.
}
\usage{
plotmix.pois(x, theta, hist = 1, comp = TRUE, h = 1, main = "",
             xlab = "Observations", ylab = "")
}
\arguments{
  \item{x}{
data, can be either a vector or a matrix with the 1st column being the observed values and the 2nd column being the corresponding frequencies.
}
  \item{theta}{
parameter values, output of emtest.pois or pmle.pois, or a vector of parameter values, include mixing proportions and component parameters.
}
  \item{hist}{
style of histogram. hist = 0,1,2 are permitted. hist = 0 means no histogram, hist = 1,2 give two styles, default value: hist = 1.
}
  \item{comp}{
a parameter specifies the form of the component probability functions in the graph.
The component probability functions will be plotted only if the order is 3 or lower.
comp = T switches on the component probability functions in the graph. If comp = F, no component probability functions will be included. 
default value: comp = T.
}
  \item{h}{
a number between 0 and 1, specifying the height of the probability function, default value: h=1 (complete probability function).
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
emtest.pois, pmle.pois, rmix.pois
}
\examples{
#generate a random sample from a 2 component Poisson mixture, 
#plot the histogram of the observations and the probability mass function with a set of 
#given parameters.
x <- rmix.pois(200,c(0.3,0.7),c(2,5))
theta <- c(0.25,0.75,1.8,5.2)
plotmix.pois(x,theta)
}
\keyword{dplot}
