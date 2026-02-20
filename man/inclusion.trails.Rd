\name{inclusion.trails}
\alias{inclusion.trails}
\docType{data}
\title{Porphyroblast inclusion trails dataset}
\description{
This dataset corresponds to two-dimensional axial orientations of porphyroblast inclusion trails, obtained from oriented rock samples of the Cenozoic Betic-Rif orogen of southern Spain and northern Morocco. About 60\% of this data was published by Aerden et al (2022) and Ruiz-Fuentes and Aerden (2023) and the rest collected subsequently by the same authors. The whole dataset was analyzed in Alonso-Pena et al (2026).
}
\usage{data(inclusion.trails)}
\format{
  A data frame with 2827 observations on the following 6 variables:
  \describe{
    \item{Tectonic.UNIT}{Factor indicating the tectonic unit to which the observation belongs.}
    \item{Area}{Factor indicating the broader geographical area.}
    \item{Subarea}{Factor specifying a subdivision within the area.}
    \item{Strike.of.Si}{Numeric variable giving the axial orientation (strike) of the inclusion trails, measured in degrees. Directions are recorded clockwise from North (0 degrees).}
    \item{Mineral}{Factor indicating the mineral type of the porphyroblast.}
   \item{Type}{Factor indicating the method used to obtain the strike measurement.}
  }
}

\details{
The dataset contains measurements of porphyroblast inclusion-trail orientations collected in a structural geology context. The key variable is \code{Strike.of.Si}, which represents axial directional data expressed in degrees. 

Because the data are axial, orientations differing by 180 degrees represent the same structural direction. All angles are measured clockwise, with 0 degrees corresponding to North.
}
\source{
Ruiz-Fuentes, A., and Aerden, D.G.A.M. (2023).
Deciphering western Mediterranean kinematics using metamorphic porphyroblasts from the Alpujarride Complex (Betic Cordillera).
\emph{Journal of Structural Geology}, 168, 104823.
\doi{10.1016/j.jsg.2023.104823}

Aerden, D.G.A.M., Farrell, T.P., Baxter, E.F., Stewart, E.M., Ruiz-Fuentes, A., and Bouybaouene, M. (2022).
Refined tectonic evolution of the Betic-Rif orogen through integrated 3-D microstructural analysis and Sm-Nd dating of garnet porphyroblasts.
\emph{Tectonics}, 41(10), e2022TC007366.
\doi{10.1029/2022TC007366}

Alonso-Pena, M., Aerden, D.G.A.M. and Angulo, J.M. (2026). Quantifying tectonic relevance of porphyroblast inclusion trails: axial mode estimation and clustering.
}
\examples{
data(inclusion.trails)
directions <- circular(inclusion.trails$Strike.of.Si,units="degrees",
                        template="geographics",modulo="pi")
}
\keyword{datasets}




