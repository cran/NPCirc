\name{pm10}
\alias{pm10}
\docType{data}
\title{
Pm10 particles in Pontevedra, Spain
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset contains measurements on pm10 particles concentration in the city of Pontevedra, Spain, as well as wind speed and wind direction recordings. The measurements were taken in a meteorological station in Pontevedra, during the year 2019, every six hours. 
}
\usage{data("pm10")}
\format{
  A data frame with 1156 observations on the following 3 variables.
  \describe{
    \item{\code{pm10}}{A numeric vector containing the concentration of pm10 particles (in micrograms per cubic meter).}
    \item{\code{direction}}{A numeric vector containing the wind direction in degrees between 0 and 360. }
    \item{\code{speed}}{A numeric vector containing the wind speed measured in kilometers per hour.}
  }
}

\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Meteogalicia \url{https://www.meteogalicia.gal/}.
}

\examples{
data(pm10)
}
\keyword{datasets}
