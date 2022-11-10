\name{dpreg.circ}
\alias{dpreg.circ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint parametric estimation of mean and dispersion functions in circular double Poisson models
}
\description{
Function \code{dpreg.circ} implements the parametric joint estimator of the mean and dispersion functions when the covariate is circular and the conditional distribution is a double Poisson, a particular case of the double exponential family. It is assumed that the logarithm of the mean and the logarithm of the dispersion are sums of sine and cosine terms.
}
\usage{
dpreg.circ(x, y, k = 2, ktilde = 1, startvmu = NULL, startvgam = NULL,
          tol= 0.000001, maxit = 300)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of data for the independent variable. The object is coerced to class \code{\link{circular}}.}
  \item{y}{Vector of data for the dependent variable. This must be same length as \code{x} and should contain counts.}
  \item{k}{Number of components for modeling the logarithm of the mean, including the intercept. Equivalent to the number of parameters to be estimated for the mean function.}
  \item{ktilde}{Number of components for modeling the logarithm of the dispersion, including the intercept. Equivalent to the number of parameters to be estimated for the dispersion function.}
  \item{startvmu}{Vector of length \code{k} containing the initial values for the parameters corresponding to the estimation of the mean.}
  \item{startvgam}{Vector of length \code{ktilde} containing the initial values for the parameters corresponding to the estimation of the dispersion.}
  \item{tol}{Tolerance parameter for convergence in the numerical estimation.}
  \item{maxit}{Maximum number of iterations in the numerical estimation.}

}
\details{
See Alonso-Pena et al. (2022) for details.
}
\value{
A list containing the following components:

\item{datax, datay}{ Original dataset.}
\item{coefficients_mu}{ A vector of length \code{k} containing the estimators for the parameters corresponding to the mean.}
\item{coefficients_mu}{ A vector of length \code{ktilde} containing the estimators for the parameters corresponding to the dispersion.}
\item{numit}{ Number of iterations needed for convergence.}
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


\examples{
\donttest{
data(spikes)
direction<-circular(spikes$direction,units="degrees")
counts<-spikes$counts
output<-dpreg.circ(direction, counts,  k = 5, ktilde = 3)

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
