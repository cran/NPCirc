\name{dragonfly}
\alias{dragonfly}
\docType{data}
\title{ Orientations of dragonflies
}
\description{
The data, presented in Batschelet (1981), consists on the orientation of 214 dragonflies with respect to the sun's azimuth.
}
\usage{data(dragonfly)}
\format{
  A single-column data frame with 214 observations in radians.
}
\details{
Data were originally recorded in degrees.
}
\source{
Batschelet, E. (1981) \emph{Circular Statistics in Biology}. Academic Press, New York.
}
\examples{
data(dragonfly)
plot(circular(dragonfly), shrink=1.3)
t <- seq(0,2*pi,length=500)
dens <- kern.den.circ(dragonfly$orientation, t, nu=10)
lines(circular(t), dens)
}
\keyword{datasets}

