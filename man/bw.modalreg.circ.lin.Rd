\name{bw.modalreg.circ.lin}
\alias{bw.modalreg.circ.lin}
\alias{bw.modalreg.lin.circ}
\alias{bw.modalreg.circ.circ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smoothing parameter selection for circular multimdodal regression
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function \code{bw.modalreg.circ.lin} provides the modal cross-validation smoothing parameters for the multimodal regression estimator when the covariate is circular and the response variable is linear.

Function \code{bw.modalreg.circ.circ} provides the modal cross-validation smoothing parameters for the multimodal regression estimator when the covariate and the response variable are circular.

Function \code{bw.modalreg.lin.circ} provides the modal cross-validation smoothing parameters for the multimodal regression estimator when the covariate is linear and the response variable is circular.
}
\usage{

bw.modalreg.circ.lin(x, y, lower = NULL, upper = NULL, maxit = 500, tol = 0.00001)
bw.modalreg.circ.circ(x, y, lower = NULL, upper = NULL, maxit = 500, tol = 0.00001)
bw.modalreg.lin.circ(x, y, lower = NULL, upper = NULL, maxit = 500, tol = 0.00001)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of data for the independent variable. The object is coerced to class circular when using functions bw.modalreg.circ.lin and bw.modalreg.circ.circ.}
  \item{y}{Vector of data for the dependent variable. This must be same length as x. The object is coerced to class circular when using functions bw.modalreg.circ.circ and bw.modalreg.lin.circ.}
  \item{lower, upper}{Vectors of length two with the \code{lower} and \code{upper} boundaries of the intervals to be used in the search for the values of the smoothing parameters.}
  \item{maxit}{Maximum number of iterations in the estimation through the conditional (circular) mean shift.}
  \item{tol}{Tolerance parameter for convergence in the estimation through the conditional (circular) mean shift.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
See Alonso-Pena and Crujeiras (2022) for details.

The NAs will be automatically removed.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A vector of length two with the first component being the value of the smoothing parameter associated to the predictor variable and with the second component being the value of the smoothing parameter associated to the response variable.
}
\references{
%% ~put references to the literature/web site here ~
Alonso-Pena, M. and Crujeiras, R. M. (2022). Analizing animal escape data with circular nonparametric multimodal regression. Annals of Applied Statistics. (To appear).
}
\author{
%%  ~~who you are~~
Maria Alonso-Pena and Rosa M. Crujeiras.
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{modalreg.circ.lin}}, \code{\link{modalreg.circ.circ}}, \code{\link{modalreg.lin.circ}}
}

\examples{
\donttest{
  # Circ-lin
  set.seed(8833)
  n1<-100
  n2<-100
  gamma<-8
  sigma<-1.5
  theta1<-rcircularuniform(n1)
  theta2<-rcircularuniform(n2)
  theta<-c(theta1,theta2)
  y1<-2*sin(2*theta1)+rnorm(n1,sd=sigma)
  y2<-gamma+2*sin(2*theta2)+rnorm(n2,sd=sigma)
  y<-as.numeric(c(y1,y2))
  bw<-bw.modalreg.circ.lin(theta, y)
  
  # Lin-circ
  n1<-100
  n2<-100
  con<-8
  set.seed(8833)
  x1<-runif(n1)
  x2<-runif(n2)
  phi1<-(6*atan(2.5*x1-3)+rvonmises(n1,m=0,k=con))%%(2*pi)
  phi2<-(pi+6*atan(2.5*x2-3)+rvonmises(n2,m=0,k=con))%%(2*pi)
  x<-c(x1,x2)
  phi<-c(phi1,phi2)
  bw<-bw.modalreg.lin.circ(x, phi)
  
  # Circ-circ
  n1<-100
  n2<-100
  con<-10
  set.seed(8833)
  theta1<-rcircularuniform(n1)
  theta2<-rcircularuniform(n2)
  phi1<-(2*cos(theta1)+rvonmises(n1,m=0,k=con))%%(2*pi)
  phi2<-(3*pi/4+2*cos(theta2)+rvonmises(n2,m=0,k=con))%%(2*pi)
  theta=c(theta1,theta2)
  phi=c(phi1,phi2)
  bw<-bw.modalreg.lin.circ(theta, phi)
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
