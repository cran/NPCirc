\encoding{latin1}
\name{circsizer.density}
\alias{circsizer.density}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ CircSiZer map for density }
\description{
   This function plots the CircSiZer map for circular density estimation based on circular kernel methods, as described in Oliveira et al. (2013). The 
   CircSiZer is an extension of SiZer proposed by Chaudhuri and Marron (1999) to circular data.
}
\usage{
circsizer.density(x, NU, ngrid=250, alpha=0.05, B=500, type=3, raw.data=FALSE,		
log.scale=TRUE, zero=pi/2, clockwise=TRUE, title=NULL, labels=NULL, 
label.pos=NULL, rad.pos=NULL)
}
\arguments{
  \item{x}{ Sample of angles in radians (between \eqn{0} and \eqn{2\pi}) from which the estimate is to be computed.}
  \item{NU}{ Vector of smoothing parameters. Values of \code{NU} must be positive. \code{NU} will be coerced to be equally spaced. Length of vector \code{NU}
	     must be at least 2.}
  \item{ngrid}{ Integer indicating the number of equally spaced angles between \eqn{0} and \eqn{2\pi} where the estimator is evaluated. Default 
                \code{ngrid=250}.}
  \item{alpha}{ Significance level for the CircSiZer map. Default \code{alpha=0.05}.}
  \item{B}{ Integer indicating the number of bootstrap samples to estimate the standard deviation of the derivative estimator. Default \code{B=500}.}
  \item{type}{ Number indicating the labels to display in the plot: 1 (directions), 2 (hours), 3 (angles in radians), 4 (angles in degrees) or 5 (months). 
               Default \code{type=3}.}
  \item{raw.data}{ Logical, if \code{TRUE}, points indicated by \code{x} are stacked on the perimeter of the circle. Default is \code{FALSE}.}
  \item{log.scale}{ Logical, if \code{TRUE}, the CircSiZer map is plotted in the scale \eqn{-log_{10}(NU)}. Default is \code{TRUE}. See Details.}
  \item{zero}{ Where to place the starting (zero) point. Defaults to the North position.}
  \item{clockwise}{ Whether to interpret positive positions as clockwise from the starting point. The default is clockwise (\code{clockwise=TRUE}).}
  \item{title}{ Title for the plot.}
  \item{labels}{ Character or expression vector of labels to be placed at the \code{label.pos}. \code{label.pos} must also be supplied.}
  \item{label.pos}{ Vector indicating the position (between \eqn{0} and \eqn{2\pi}) at which the labels are to be drawn.}
  \item{rad.pos}{ Vector (between \eqn{0} and \eqn{2\pi}) with the drawing position for the radius.}
}
\details{
With CircSiZer, significance features (peaks and valleys) in the data are seeked via the construction of confidence intervals for the scale-space version of 
the smoothed derivative curve, as it is described in Oliveira et al. (2013). Thus, for a given point and a given value of the smoothing parameter, the curve 
is significantly increasing (decreasing) if the confidence interval is above (below) 0 and if the confidence interval contains 0, the curve for that value
of the smoothing parameter and at that point does not have a statistically significant slope. This information is displayed in a circular color map, the 
CircSiZer map, in such a way that, at a given point, the performance of the estimated curve is represented by a color ring with radius proportional to the 
value of the smoothing parameter.

Differents colors allow to  indentify peaks and valleys. Blue color indicates locations where the curve is significantly increasing; red color shows where 
it is significantly decreasing and purple indicates where it is not significantly different from zero. Gray color correspons to those regions where there is 
not enough data to make statements about significance. Thus, at a given bandwidth, a significant peak can be identified when a region of significant 
positive gradient is followed by a region of significant negative gradient (i.e. blue-red pattern), and a significant trough by the reverse (red-blue 
pattern), taking clockwise as the positive sense of rotation.

If \code{log.scale=TRUE} then, the values of the considered smoothing parameters \code{NU} are transformed to \eqn{-log_{10}} scale, i.e, a sequence of equally
spaced smoothing parameters according to the parameters 			
\code{-log10(max(NU))}, \code{-log10(min(NU))} and \code{length(NU)} is used. Hence, small values 
of this parameter corresponds with larger rings and large values corresponds with smaller rings. Whereas if \code{log.scale=FALSE}, small values of this 
parameter corresponds with smaller rings and large values corresponds with larger rings.

The NAs will be automatically removed.
}
\value{
  CircSiZer map for density.
}
\references{ 
Chaudhuri, P. and Marron, J.S. (1999). SiZer for exploration of structures in curves, \emph{Journal of the American Statistical Association}, \bold{94}, 
807--823.

Oliveira, M., Crujeiras, R.M. and Rodríguez--Casal (2013) CircSiZer: an exploratory tool for circular data. http://arxiv.org/abs/1210.2555v2.
}


\author{ María Oliveira, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\examples{
# set.seed(2012)
# x <- rcircmix(100,model=7)
# circsizer.density(x, NU=seq(0,50,length=12),type=4,zero=0,clockwise=FALSE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ circular density }% __ONLY ONE__ keyword per line
