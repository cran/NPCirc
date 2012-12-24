\name{speed.wind}
\alias{speed.wind}
\alias{speed.wind2}
\docType{data}
\title{ Wind speed and wind direction data
}
\description{
This dataset consists of hourly observations of wind direction and wind speed in winter season (from November to February)
from 2003 until 2012 in the Atlantic coast of Galicia (NW--Spain). Data are registered by a buoy located at longitude -0.210E and latitude 43.500N in the 
Atlantic Ocean. The dataset \code{speed.wind2}, analized in Oliveira et al. (2013), is a subset of \code{speed.wind} which is obtained by taking the 
observations with a lag period of 95 hours.
}
\usage{
data(speed.wind)
data(speed.wind2)
}
\format{
  \code{speed.wind} is a data frame with 19488 observations on six variables: day, month, year, hour, wind speed (in m/s) and wind direction 
  (in degrees). \code{speed.wind2} is a subset with 200 observations.
}
\details{
Data contains NAs.
There is no data in November 2005, December 2005, January 2006, February 2006, February 2007, February 2009 and November 2009. Months of November 2004, 
December 2004, January 2007, December 2009 are not fully observed.
}
\source{
Data can be freely downloaded from the Spanish Portuary Authority (http://www.puertos.es).
}
\references{ Oliveira, M., Crujeiras, R.M. and Rodríguez--Casal (2013) CircSiZer: an exploratory tool for circular data. http://arxiv.org/abs/1210.2555v2.
}
\examples{
data(speed.wind2)

t <- seq(0,2*pi,length=500)
dir <- rad(speed.wind2$Direction)

# Density
plot(circular(dir), stack=TRUE, zero=pi/2, rotation="clock", axes=FALSE)
axis.circular(at=circular(seq(0,2*pi,by=pi/4)),
labels=c("N","NE","E","SE","S","SW","W","NW","N"), zero=pi/2, rotation="clock")
rose.diag(circular(dir), bins=16, add=TRUE, axes=FALSE, zero=pi/2, rotation="clock")

rose.diag(circular(dir), bins=16, shrink=1.1, axes=FALSE, zero=pi/2, rotation="clock")
axis.circular(at=circular(seq(0,7*pi/4,by=pi/4)),
labels=c("N","NE","E","SE","S","SW","W","NW"), zero=pi/2, rotation="clock")
lines(circular(t), kern.den.circ(dir,t,nu=1), lwd=2, lty=2, zero=pi/2, rotation="clock")
lines(circular(t), kern.den.circ(dir,t,nu=10), lwd=2, lty=1, zero=pi/2, rotation="clock")
lines(circular(t), kern.den.circ(dir,t,nu=60), lwd=2, lty=3, zero=pi/2, rotation="clock")

# Regression
vel <- speed.wind2$Speed
nas <- which(is.na(vel))
dir <- dir[-nas]
vel <- vel[-nas]
radial.plot(vel, dir, rp.type="s", start=pi/2, clockwise=TRUE, point.col="gray",
radial.lim=c(0,15), label.pos=seq(0,7*pi/4,by=pi/4),
labels=c("N","NE","E","SE","S","SW","W","NW"))
radial.plot(as.numeric(kern.reg.circ(dir,vel,t,nu=1,method="LL")), t, rp.type="p",
add=TRUE, start=pi/2, clockwise=TRUE, radial.lim=c(0,15), lwd=2, lty=2)
radial.plot(as.numeric(kern.reg.circ(dir,vel,t,nu=10,method="LL")),t, rp.type="p",
add=TRUE, start=pi/2, clockwise=TRUE, radial.lim=c(0,15), lwd=2, lty=1)
radial.plot(as.numeric(kern.reg.circ(dir,vel,t,nu=60,method="LL")),t, rp.type="p",
add=TRUE, start=pi/2, clockwise=TRUE, radial.lim=c(0,15), lwd=2, lty=3)

}
\keyword{datasets}

