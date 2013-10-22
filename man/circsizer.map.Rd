\encoding{latin1}
\name{circsizer.map}
\alias{circsizer.map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ CircSiZer map }
\description{
   This function plots the CircSiZer map for \code{circsizer} objects.
}
\usage{
circsizer.map(circsizer.object, type, zero, clockwise, title=NULL, labels=NULL, 
label.pos=NULL, rad.pos=NULL, raw.data=FALSE)
}

\arguments{
  \item{circsizer.object}{ An object of class \code{circsizer}, i.e., output from functions \code{circsizer.density} or \code{circsizer.regression}.}
  \item{type}{ Number indicating the labels to display in the plot: 1 (directions), 2 (hours), 3 (angles in radians), 4 (angles in degrees) or 5 (months).}
  \item{zero}{ Where to place the starting (zero) point.}
  \item{clockwise}{ Whether to interpret positive positions as clockwise from the starting point.}
  \item{title}{ Title for the plot.}
  \item{labels}{ Character or expression vector of labels to be placed at the \code{label.pos}. \code{label.pos} must also be supplied.}
  \item{label.pos}{ Vector indicating the position (between \eqn{0} and \eqn{2\pi}) at which the labels are to be drawn.}
  \item{rad.pos}{ Vector (between \eqn{0} and \eqn{2\pi}) with the drawing position for the radius.}
  \item{raw.data}{ Logical, if \code{TRUE}, points indicated by \code{x} are stacked on the perimeter of the circle. Default is \code{FALSE}.}
}


\value{
  CircSiZer map.
}

%\references{ 
%Chaudhuri, P. and Marron, J.S. (1999). SiZer for exploration of structures in curves, \emph{Journal of the American Statistical Association}, \bold{94}, 
%807--823.
%
%Oliveira, M., Crujeiras, R.M. and Rodríguez--Casal (2013) CircSiZer: an exploratory tool for circular data. \emph{Environmental and Ecological Statistics}, 
%DOI: 10.1007/s10651-013-0249-0.
%}

\author{ María Oliveira, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\seealso{\code{\link{circsizer.density}}, \code{\link{circsizer.regression}} }


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ circsizer }% __ONLY ONE__ keyword per line
