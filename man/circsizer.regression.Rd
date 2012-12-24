\encoding{latin1}
\name{circsizer.regression}
\alias{circsizer.regression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ CircSiZer map for regression }
\description{
   This function plots the CircSiZer map for circular regression estimation based on circular kernel methods, as described in Oliveira et al. (2013). The 
   CircSiZer is an extension of SiZer proposed by Chaudhuri and Marron (1999) to circular data.
}
\usage{
circsizer.regression(x, y, NU, ngrid=150, alpha=0.05, B=500, B2=250, type=3, 
log.scale=TRUE, zero=pi/2, clockwise=TRUE, title=NULL, labels=NULL, 
label.pos=NULL, rad.pos=NULL)
}
\arguments{
  \item{x}{ Sample of angles in radians (between \eqn{0} and \eqn{2\pi}) for the circular covariate.}
  \item{y}{ Sample of angles for the dependent linear variable. This must be same length as \code{x}.}
  \item{NU}{ Vector of smoothing parameters. Values of \code{NU} must be positive. \code{NU} will be coerced to be equally spaced. Length of vector \code{NU}
	     must be at least 2.}
  \item{ngrid}{ Integer indicating the number of equally spaced angles between \eqn{0} and \eqn{2\pi} where the estimator is evaluated. Default 
                \code{ngrid=150}.}
  \item{alpha}{ Significance level for the CircSiZer map. Default \code{alpha=0.05}.}
  \item{B}{ Integer indicating the number of bootstrap samples to estimate the standard deviation of the derivative estimator. Default \code{B=500}.}
  \item{B2}{ Integer indicating the number of bootstrap samples to compute the denominator in Step 2 of algorithm described in Oliveira et al. (2013). 
             Default \code{B=250}.}
  \item{type}{ Number indicating the labels to display in the plot: 1 (directions), 2 (hours), 3 (angles in radians), 4 (angles in degrees) or 5 (months). 
               Default \code{type=3}.}
  \item{log.scale}{ Logical, if \code{TRUE}, the CircSiZer map is plotted in the scale \eqn{-log_{10}(NU)}. Default is \code{TRUE}.}
  \item{zero}{ Where to place the starting (zero) point. Defaults to the North position.}
  \item{clockwise}{ Whether to interpret positive positions as clockwise from the starting point. The default is clockwise (\code{clockwise=TRUE}).}
  \item{title}{ Title for the plot.}
  \item{labels}{ Character or expression vector of labels to be placed at the \code{label.pos}. \code{label.pos} must also be supplied.}
  \item{label.pos}{ Vector indicating the position (between \eqn{0} and \eqn{2\pi}) at which the labels are to be drawn.}
  \item{rad.pos}{ Vector (between \eqn{0} and \eqn{2\pi}) with the drawing position for the radius.}
}
\details{
  See Details Section of \code{\link{circsizer.density}}.
  The NAs will be automatically removed.
}
\value{
  CircSiZer map for regression.
}
\references{
Chaudhuri, P. and Marron, J.S. (1999). SiZer for exploration of structures in curves, \emph{Journal of the American Statistical Association}, \bold{94}, 
807--823.

Oliveira, M., Crujeiras, R.M. and Rodríguez--Casal (2013) CircSiZer: an exploratory tool for circular data. http://arxiv.org/abs/1210.2555v2.
}


\author{ María Oliveira, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\examples{
# Not run: the code works but it is slow
# set.seed(2012)
# n <- 100
# x <- seq(0,2*pi,length=n)
# y <- sin(x)+sqrt(0.5)*rnorm(n)
# circsizer.regression(x, y, NU=seq(10,60,by=5), title="CircSiZer for regression")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ circular regression }% __ONLY ONE__ keyword per line
