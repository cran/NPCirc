\name{modalreg.circ.lin}
\alias{modalreg.circ.lin}
\alias{modalreg.circ.circ}
\alias{modalreg.lin.circ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Circular multimodal regression estimation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function \code{modalreg.circ.lin} implements the nonparametric multimodal regression estimator for a circular covariate and a real-valued response, as described in Alonso-Pena and Crujeiras (2022). It takes the von Mises distribution as the kernel associated to the predictor variable and the normal distribution as the kernel associated to the response variable.

Function \code{modalreg.circ.circ} implements the nonparametric multimodal regression estimator for a circular covariate and a circular response, as described in Alonso-Pena and Crujeiras (2022). It takes the von Mises distribution as the kernel associated to the predictor variable and the response variables.

Function \code{modalreg.lin.circ} implements the nonparametric multimodal regression estimator for a real-valued covariate and a circular response, as described in Alonso-Pena and Crujeiras (2022). It takes the normal distribution as the kernel associated to the predictor variable and the von Mises distribution as the kernel associated to the response variable.

}

\usage{
modalreg.circ.lin(x, y, t=NULL, bw=NULL, tol = 0.0001, maxit = 500,
    from = circular(0),to = circular(2 * pi), len = 300)

modalreg.circ.circ(x, y, t=NULL, bw=NULL, tol = 0.00001,
    maxit = 500, from = circular(0), to = circular(2 * pi), len = 300)

modalreg.lin.circ(x, y, t=NULL, bw=NULL, tol = 0.0001, maxit = 500, len=300)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of data for the independent variable. The object is coerced to class \code{\link{circular}} when using functions \code{modalreg.circ.lin} and \code{modalreg.circ.circ}.}
  \item{y}{Vector of data for the dependent variable. This must be same length as \code{x}. The object is coerced to class \code{\link{circular}} when using functions \code{modalreg.circ.circ} and \code{modalreg.lin.circ}.}
  \item{t}{Points where the regression function  is estimated. If \code{NULL}, equally spaced points are used according to the parameters \code{from}, \code{to} and \code{len}.}
  \item{bw}{Vector of length two with the values of the smoothing parameters to be used. The first component corresponds to the smoothing parameter associated to the predictor variable and the second component is the parameter associated to the response variable. If \code{NULL}, the parameters are selected via modal cross-validation.}
  \item{tol}{Tolerance parameter for convergence in the estimation through the conditional (circular) mean shift.}
  \item{maxit}{Maximum number of iterations in the estimation through the conditional (circular) mean shift.}
  \item{from, to}{ Left and right-most points of the grid at which the regression multifunction is to be estimated. The objects are coerced to class \code{\link{circular}}.}
  \item{len}{ Number of equally spaced points at which the regression multifunction is to be estimated.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
See Alonso-Pena and Crujeiras (2022) for details.

The NAs will be automatically removed.
}
\value{

A list containing the following components:

\item{datax, datay}{ Original dataset.}
\item{x}{ The n coordinates of the points where the regression multifunction is estimated.}
\item{y}{ A list with dimension the length of the number of evaluation points  containing  the estimated values of the multidunfunction for each evaluation point.}
\item{bw}{ A verctor of lenght two with the smoothing parameters used.}
\item{n}{ The sample size after elimination of missing values.}
\item{call}{ The call which produced the result.}
\item{data.name}{ The deparsed name of the x argument.}
\item{has.na}{ Logical, for compatibility (always FALSE).}
}
\references{
%% ~put references to the literature/web site here ~
Alonso-Pena, M. and Crujeiras, R. M. (2022). Analizing animal escape data with circular nonparametric multimodal regression. Annals of Applied Statistics. (To appear).
}
\author{
%%  ~~who you are~~
Maria Alonso-Pena and Rosa M. Crujeiras.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{bw.modalreg.circ.lin}}, \code{\link{bw.modalreg.circ.circ}}, \code{\link{bw.modalreg.lin.circ}}
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
  fit<-modalreg.circ.lin(theta,y,bw=c(10,1.3))
  
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
  fit<-modalreg.lin.circ(x, phi, bw=c(0.1,2.5))
  
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
  fit<-modalreg.circ.circ(theta, phi, bw=c(30,3))
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
