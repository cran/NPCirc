\name{flywheels}
\alias{flywheels}
\docType{data}
\title{
Flywheel measurements
}
\description{
The data consists on measurements on mechanical flywheels, which are tools to regulate an engine's rotation in the automotive industry. The dataset contains 60 observations, containing the flywheels' angle of imbalance and the magnitude required to correct them. Four different types of metals were employed in the production process. 
}
\usage{data("flywheels")}
\format{
  A data frame with 60 observations on the following 3 variables.
  \describe{
    \item{\code{angle}}{A numeric vector containing the angles of imbalance (in radians)}
    \item{\code{weight}}{A numeric vector containing the magnitudes of correction (in inch-ounces).}
    \item{\code{group}}{A vector with numbers from 1 to 4 indicating the type of metal.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The data were analyzed in Anderson-Cook (1999) and Alonso-Pena et al. (2021).
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
The dataset was obtained from Anderson-Cook (1999).
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Anderson-Cook, C.M. (1999). A tutorial on one-way analysis of circular-linear data. \emph{Journal of Quality Technology}, \strong{31}(1), 109--119.

Alonso-Pena, M., Ameijeiras-Alonso, J. and Crujeiras, R.M. (2021) Nonparametric tests for circular regression. \emph{Journal of Statistical Computation and Simulation}, \strong{91}(3), 477--500.
}
\examples{
data(flywheels)
}
\keyword{datasets}
