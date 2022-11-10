\name{spikes}
\alias{spikes}
\docType{data}
\title{
Neuronal spikes in a macaque monkey
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset, given by Kohn and Movshon (2003), contains observations from an experimental study on macaque monkeys, where a macaque received a visual stimulus. The dataset is composed of 68 observations, consisting on the angle of stimulus and the number of spikes produced in a MT/V5 neuron.
}
\usage{data("spikes")}
\format{
  A data frame with 68 observations on the following 2 variables.
  \describe{
    \item{\code{direction}}{A numeric vector containing the direction of the visual stimulus (in degrees).}
    \item{\code{counts}}{A numeric vector containing the number of spikes.}
  }
}
\details{
See Kohn and Movshon (2003) for details. This dataset was also studied in Alonso-Pena et al (2022). 
}
\source{
The dataset was obtained by Kohn and Movshon (2003). 
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Kohn, A., and Movshon, J.A. (2003). Neuronal adaptation to visual motion in area {MT} of the macaque. \emph{Neuron}, \strong{39}(4), 681--691.

Alonso-Pena, M., Gijbels, I. and Crujeiras, R.M. (2022). Flexible joint modeling of mean and dispersion for the directional tuning of neuronal spike counts. \emph{Under review}.
}
\examples{
data(spikes)
}
\keyword{datasets}
