\name{temp.wind}
\alias{temp.wind}
\docType{data}
\title{ Temperature and wind direction data
}
\description{
The data, analized in Oliveira et al. (2013), consists of observations of temperature and wind direction during the austral summer season 2008-2009 
(from November 2008 to March 2009) in Vinciguerra (Tierra del Fuego, Argentina).
}
\usage{data(temp.wind)}
\format{
  A data frame with 3648 observations on four variables: \code{Date}, \code{Time}, \code{Temperature} (in degrees Celsius) and \code{Direction} (in degrees).
}
\details{
Data contains NAs.
}
\source{
The authors want to acknowledge Prof. Augusto Pérez--Alberti for providing the data, collected within the Project POL2006-09071 from the Spanish Ministry 
of Education and Science.
}
\references{ Oliveira, M., Crujeiras R.M. and Rodríguez--Casal, A. (2013) Nonparametric circular methods for exploring environmental data. 
\emph{Environmental and Ecological Statistics}, \bold{20}, 1--17.
}
\examples{
data(temp.wind)

winddir <- rad(temp.wind[,4]) # wind direction in radians
temp <- temp.wind[,3]
nas <- which(is.na(winddir))
winddir <- winddir[-nas] 
temp <- temp[-nas]

# value of the smoothing parameter selected using the function nu.LSCV.reg 
# with method="LL"
nu <- 3.41 
t <- seq(0,2*pi,length=100)
est <- kern.reg.circ(winddir, temp, t, nu=nu, method="LL")

# Circular representation
radial.plot(temp, winddir, rp.type="s", labels=c("N","NE","E","SE","S","SW","W","NW"),
start=pi/2, clockwise=TRUE, label.pos=seq(0,7*pi/4,by=pi/4), lwd=2, point.col="grey50",
radial.lim=c(-10,15))
radial.plot(as.vector(est), rp.type="p", start=pi/2, clockwise=TRUE, lwd=2,
radial.lim=c(-10,15), add=TRUE)

# Linear representation
plot(t, est, type="l", xlab="", ylab="Temperatute (ºC)", axes=FALSE)
axis(1, labels=c("N","NE","E","SE","S","SW","W","NW","N"), at=seq(0,2*pi,by=pi/4))
axis(2)
}
\keyword{datasets}

