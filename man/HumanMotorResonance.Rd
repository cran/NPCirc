\name{HumanMotorResonance}
\alias{HumanMotorResonance}
\docType{data}
\title{
Human motor resonance data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset, given by Puglisi et al. (2017), contains measurements of motor resonance responses in humans, who were requested to observe a movement of a rythmic hand flexion-extension. The dataset is composed of 70 observations, including the angular position of the hand and the resonance response.
}
\usage{data("HumanMotorResonance")}
\format{
  A data frame with 70 observations on the following 2 variables.
  \describe{
    \item{\code{Angular.position}}{A numeric vector containing the angular position of the hand (in radians).}
    \item{\code{Reflex.Amplitude}}{A numeric vector containing the resonance response.}}
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
See Puglisi et al. (2017) for details.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
The dataset was obtained by Puglisi et al. (2017).
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Puglisi, G.,  Leonetti, A., Landau, A., Fornia, L., Cerri, G. and Borroni, P. (2017). The role of attention in human motor resonance. \emph{PLOS ONE}, \strong{12}(5), e0177457.
}
\examples{
data(HumanMotorResonance)
}
\keyword{datasets}
