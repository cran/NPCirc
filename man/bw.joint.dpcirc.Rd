\name{bw.joint.dpcirc}
\alias{bw.joint.dpcirc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smoothing parameter selection for circular double Poisson regression
}
\description{
Function \code{bw.joint.dpcirc} provides the smoothing parameters for the nonparametric joint estimator of the mean and dispersion functions when the conditional density is a double Poisson. It performs a joint cross-validation search. 
}
\usage{
bw.joint.dpcirc(x, y, startvmu = NULL, startvgam = NULL, lower=c(0.05,0.05), 
  upper=c(50,7),tol = 0.00001, maxit = 300)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{Vector of data for the independent variable. The object is coerced to class circular.}
  \item{y}{Vector of data for the dependent variable. This must be same length as x and should contain counts.}
  \item{startvmu}{Vector of length two containing the initial values for the parameters corresponding to the estimation of the mean.}
  \item{startvgam}{Vector of length two containing the initial values for the parameters corresponding to the estimation of the dispersion.}
  \item{lower, upper}{Vectors of length two with the \code{lower} and \code{upper} boundaries of the intervals to be used in the search for the values of the smoothing parameters. The first component of each corresponds to the parameter associated to the estimation of the mean, while the second component corresponds to the estimation of the dispersion.}
  \item{tol}{Tolerance parameter for convergence in the numerical estimation.}
  \item{maxit}{Maximum number of iterations in the numerical estimation.}
}
\details{
See Alonso-Pena et al. (2022) for details.

The NAs will be automatically removed.
}
\value{
A vector of length two with the first component being the value of the smoothing parameter associated to the mean estimation and with the second component being the value of the smoothing parameter associated to the dispersion estimation.
}
\references{
Alonso-Pena, M., Gijbels, I. and Crujeiras, R.M. (2022). Flexible joint modeling of mean and dispersion for the directional tuning of neuronal spike counts. \emph{Under review}.
}
\author{
Maria Alonso-Pena, Irene Gijbels and Rosa M. Crujeiras
}

\seealso{
\code{\link{kern.dpreg.circ}}
}
\examples{
\donttest{
data(spikes)
direction<-circular(spikes$direction,units="degrees")
counts<-spikes$counts
bw.joint.dpcirc(direction, counts, lower=c(0.5,0.5), upper=c(50,7))
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
