\name{zebrafish}
\alias{zebrafish}
\docType{data}
\title{
Zebrafish
}
\description{
The data consists on measurements from an experimental study on larval zebrafish, which were startled by an imitating predator consisting on a robot disguised as an adult zebrafish. The dataset includes 502 observations corresponding to the escape directions of each fish and the angles in which the zebrafish perceive the threat.
}
\usage{data("zebrafish")}
\format{
  A data frame with 502 observations on the following 2 variables.
  \describe{
    \item{\code{stimulus}}{A numeric vector containing the angles of stimulus (in radians)}
    \item{\code{res_angle}}{A numeric vector containing the directions of escape (in radians)}
  }
}
\details{
Analysis of the zebrafish data with circular regression methods can be seen in Alonso-Pena et al. (2022).
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
The data were obtained from the Dryad Digital Repository \doi{https://doi.org/10.5061/dryad.47mq9} and first analyzed in Nair et al. (2017).
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Nair, A., Changsing, K., Stewart, W.J. and McHenry, M.J. (2017). Fish prey change strategy with the direction of a threat. \emph{Proceedings of the Royal Society B: Biological Sciences}, \bold{284}, 20170393.

Alonso-Pena, M. and Crujeiras, R. M. (2022). Analizing animal escape data with circular nonparametric multimodal regression. Annals of Applied Statistics. (To appear).
}
\examples{
data(zebrafish)
}
\keyword{datasets}
