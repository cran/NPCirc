\name{bw.circ.local.lik}
\alias{bw.circ.local.lik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smoothing parameter selection for circular local likelihood regression
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function \code{bw.circ.local.lik} computes values of the smoothing (concentration) parameter for local likelihood regression when the predictor is circular and the conditional density is either gaussian, Bernoulli, Poisson or gamma. The smoothing parameter can be selected by the refined rule, the CRSC rule or cross-validation, as described in Alonso-Pena et al. (2022).
}
\usage{
bw.circ.local.lik(x, y, t, rule = NULL, p, family, startv = NULL,
  lower = 0, upper = 50, lower_ast = 0, upper_ast = 15, tol = 0.00001,
  maxit = 300, from = circular(0), to = circular(2 * pi), len = 250)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of data for the independent variable. The object is coerced to class \code{\link{circular}}.}
  \item{y}{Vector of data for the dependent variable. This must be same length as \code{x}.}
  \item{t}{ Points where the regression function  is estimated. If \code{NULL}, equally spaced points are used according to the parameters \code{from}, \code{to} and \code{len}.}
  \item{rule}{Character string giving the rule to be used to select the smoothing (concentration) parameter. This must be one of \code{"refined"} (only for \code{p = 1}), \code{"CRSC"} or \code{"cv"}.}
  \item{p}{Degree of the local sine-polynomial to be used in the estimation process. It must be 1 or 3.}
  \item{family}{Character string indicating the conditional density to be used. It must be one of \code{"gaussian"}, \code{"bernoulli"}, \code{"poisson"} or \code{"gamma"}.}
  \item{startv}{Vector containing the initial values for the estimation algorithm if family is set as \code{"bernoulli"}, \code{"poisson"} or \code{"gamma"}. The vector must be of length 2 if \code{p = 1} and length 4 if \code{p = 3}. If \code{NULL}, the initial parameters are the ones corresponding to the global mean of the responses.}
  \item{lower, upper}{Lower and upper boundary of the interval to be used in the search for the value of the smoothing parameter. Default \code{lower = 0} and \code{upper = 50}.}
  \item{lower_ast, upper_ast}{Lower and upper boundary of the interval to be used in the search for the value of the pilot smoothing parameter in the refined rule. Default \code{lower_ast = 0} and \code{upper_ast = 15}.}
  \item{tol}{Tolerance parameter for convergence in the numerical estimation. Only needed if family is one of \code{"bernoulli"}, \code{"poisson"} or \code{"gamma"}. Default is \code{tol = 0.00001}.}
  \item{maxit}{Maximum number of iterations in the numerical estimation. Only needed if family is one of \code{"bernoulli"}, \code{"poisson"} or \code{"gamma"}. Default is \code{maxit = 300}.}
  \item{from, to}{Left and right-most points of the grid at which the density is to be estimated. The objects are coerced to class \code{\link{circular}}.}
  \item{len}{Number of equally spaced points at which the density is to be estimated.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
For the \code{refined} rule, which is only available for \code{p = 1}, first a pilot concentration parameter is selected with the (E)CRSC rule using a sine-polynomial of degree 3, where the search is conducted between the values \code{lower_ast} and \code{upper_ast}. With the pilot smoothing parameter, the MISE of the estimator is approximated and the refined rule selects the parameter which minimizes the approximated MISE.

The \code{CRSC} rule selects the parameter minimizing the Circular Residual Squares Criterion if \code{family = "gaussian"} and the Extended Circular Residual Squares Criterion in the other cases.

The \code{cv} rule performs a cross-validation search.

See Alonso-Pena et al. (2022) for more details.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Value of the smoothing parameter.
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
\code{\link{circ.local.lik}}
}
\examples{
\donttest{
data(spikes)
direction<-circular(spikes$direction,units="degrees")
counts<-spikes$counts
bw.circ.local.lik(direction, counts, rule = "refined", p=1, family="poisson")
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
