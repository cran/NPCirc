\name{circ.local.lik}
\alias{circ.local.lik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Local likelihood estimation for regression with circular covariate
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function \code{circ.local.lik} computes a nonparametric estimation of a curve of interest, regarded as a transformation of the mean regression function, when the predictor is circular and the conditional density is either gaussian, Bernoulli, Poisson or gamma. It also computes the derivatives of the function of interest. It uses the method described in Alonso-Pena et al. (2022).
}
\usage{
circ.local.lik(x, y, t = NULL, bw = NULL, family, p = 1,
  startv = NULL, tol = 0.00001, maxit = 300, from = circular(0),
  to = circular(2 * pi),len = 250)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of data for the independent variable. The object is coerced to class \code{\link{circular}}.}
  \item{y}{Vector of data for the dependent variable. This must be same length as \code{x}.}
  \item{t}{ Points where the regression function  is estimated. If \code{NULL}, equally spaced points are used according to the parameters \code{from}, \code{to} and \code{len}.}
  \item{bw}{Value of the smoothing (concentration) parameter used. The value of the smoothing parameter can be chosen by using \code{bw.circ.local.lik}.}
  \item{family}{Character string indicating the conditional density to be used. It must be one of \code{"gaussian"}, \code{"bernoulli"}, \code{"poisson"} or \code{"gamma"}. When \code{family = "gaussian"}, the conditional mean is estimated; when \code{family = "bernoulli"}, the logit function is estimated and when \code{family = "poisson"} or \code{family = "gamma"}, the function of interest is the logarithm of the conditional mean.}
  \item{p}{Degree of the local sine-polynomial to be used in the estimation process. It must be 1 or 3.}
  \item{startv}{Vector containing the initial values for the estimation algorithm if family is set as \code{"bernoulli"}, \code{"poisson"} or \code{"gamma"}. The vector must be of length 2 if \code{p = 1} and length 4 if \code{p = 3}. If \code{NULL}, the initial parameters are the ones corresponding to the global mean of the responses.}
   \item{tol}{Tolerance parameter for convergence in the numerical estimation. Only needed if family is one of \code{"bernoulli"}, \code{"poisson"} or \code{"gamma"}. Default is \code{tol = 0.00001}.}
  \item{maxit}{Maximum number of iterations in the numerical estimation. Only needed if family is one of \code{"bernoulli"}, \code{"poisson"} or \code{"gamma"}. Default is \code{maxit = 300}.}
  \item{from, to}{Left and right-most points of the grid at which the density is to be estimated. The objects are coerced to class \code{\link{circular}}.}
  \item{len}{Number of equally spaced points at which the density is to be estimated.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
See Alonso-Pena et al. (2022).

The NAs will be automatically removed.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list containing the following components:

\item{datax, datay}{ Original dataset.}
\item{x}{ The n coordinates of the points where the regression function and its derivatives are estimated.}
\item{y}{ A list containing the estimated values of the function of interest and its derivatives up to order \code{p}.}
\item{bw}{ The smoothing parameter used.}
\item{n}{ The sample size after elimination of missing values.}
\item{call}{ The call which produced the result.}
\item{data.name}{ The deparsed name of the x argument.}
\item{has.na}{ Logical, for compatibility (always FALSE).}
}
\references{
%% ~put references to the literature/web site here ~
Alonso-Pena, M., Gijbels, I. and Crujeiras, R.M. (2022). A general framework for circular local likelihood regression. \emph{Under review}.
}
\author{
%%  ~~who you are~~
Maria Alonso-Pena, Irene Gijbels and Rosa M. Crujeiras.
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{bw.circ.local.lik}}
}
\examples{
\donttest{
data(spikes)
direction<-circular(spikes$direction,units="degrees")
counts<-spikes$counts
circ.local.lik(direction, counts, bw=8, p=1, family="poisson")
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
