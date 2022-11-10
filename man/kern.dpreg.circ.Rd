\name{kern.dpreg.circ}
\alias{kern.dpreg.circ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint kernel estimation of mean and dispersion functions in circular double Poisson models
}
\description{
Function \code{kern.dpreg.circ} implements the nonparametric joint estimator of the mean and dispersion functions when the covariate is circular and the conditional distribution is a double Poisson, a particular case of the double exponential family. It takes the von Mises distribution as the kernel employed for both the estimation of the mean and the dispersion. The employed estimator is a local-linear type.
}
\usage{
kern.dpreg.circ(x, y, t = NULL, bw, startvmu = NULL, startvgam = NULL, tol= 0.000001,
  maxit = 300, from =     circular(0),to = circular(2 * pi), len = 250)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of data for the independent variable. The object is coerced to class \code{\link{circular}}.}
  \item{y}{Vector of data for the dependent variable. This must be same length as \code{x} and should contain counts.}
  \item{t}{ Points where the regression function  is estimated. If \code{NULL}, equally spaced points are used according to the parameters \code{from}, \code{to} and \code{len}.}
  \item{bw}{Numeric vector of length two containing the smoothing (concentration) parameters used. The first component is the parameter used for the estimation of the mean, while the second component is used for estimating the dispersion. The value of the smoothing parameters can be chosen by using \code{bw.joint.dpcirc}.}
   \item{startvmu}{Vector of length two containing the initial values for the parameters corresponding to the estimation of the mean.}
  \item{startvgam}{Vector of length two containing the initial values for the parameters corresponding to the estimation of the dispersion.}
  \item{tol}{Tolerance parameter for convergence in the numerical estimation.}
  \item{maxit}{Maximum number of iterations in the numerical estimation.}
  \item{from, to}{Left and right-most points of the grid at which the density is to be estimated. The objects are coerced to class \code{\link{circular}}.}
  \item{len}{Number of equally spaced points at which the function is to be estimated.}

}
\details{
See Alonso-Pena et al. (2022) for details.

The NAs will be automatically removed.
}
\value{
A list containing the following components:

\item{datax, datay}{ Original dataset.}
\item{x}{ The n coordinates of the points where the regression function and its derivatives are estimated.}
\item{estim}{ A list containing the estimated values of the logarithm of the mean function (first component) and the logarithm of the dispersion function (second component).}
\item{bw}{ The vector of smoothing parameters used.}
\item{n}{ The sample size after elimination of missing values.}
\item{call}{ The call which produced the result.}
\item{data.name}{ The deparsed name of the x argument.}
\item{has.na}{ Logical, for compatibility (always FALSE).}
}
\references{
Alonso-Pena, M., Gijbels, I. and Crujeiras, R.M. (2022). Flexible joint modeling of mean and dispersion for the directional tuning of neuronal spike counts. \emph{Under review}.
}
\author{
Maria Alonso-Pena, Irene Gijbels and Rosa M. Crujeiras
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{bw.joint.dpcirc}}
}
\examples{
\donttest{
data(spikes)
direction<-circular(spikes$direction,units="degrees")
counts<-spikes$counts
output<-kern.dpreg.circ(direction, counts, bw=c(7.41,4.47))

 # Plot the data
plot(as.numeric(output$datax),output$datay,pch=16,xlab="Stimulus direction",
      ylab="Spike number")
# Represent the estimated mean
points(as.numeric(output$x),exp(output$estim[[1]]),type="l",col=2,lwd=2)

# Represent the estimated dispersion
plot(as.numeric(output$x),exp(output$estim[[2]]),type="l",col=2,lwd=2,
      xlab="Stimulus direction",ylab="Dispersion")
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
