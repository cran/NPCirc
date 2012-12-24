\encoding{latin1}
\name{kern.reg.circ}
\alias{kern.reg.circ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Nonparametric circular kernel regression estimation }
\description{
    This function implements the Nadaraya-Watson estimator and the Local-Linear estimator for circular-linear data (circular covariate and linear response),
    as described in Di Marzio et al. (2009) and Oliveira et al. (2013), taking the von Mises distribution as kernel.
}
\usage{
kern.reg.circ(x, y, t=NULL, nu, method="LL", tol=300, from=0, to=2*pi, len=250)
}
\arguments{
  \item{x}{ Sample of angles in radians (between \eqn{0} and \eqn{2\pi}) for the circular covariate.}
  \item{y}{ Sample of angles for the dependent linear variable. This must be same length as \code{x}.}
  \item{t}{ Vector of angles (in radians) where to evaluate the estimator. If \code{NULL} equally spaced points are used according to the parameters 
            \code{from}, \code{to} and \code{len}.}
  \item{nu}{ Smoothing parameter to be used. The value of \code{nu} can be chosen by using the function \code{\link{nu.LSCV.reg}}}
  \item{method}{ Character string giving the estimator to be used. This must be one of \code{"LL"} for Local-Linear estimator or \code{"NW"} for 
                 Nadaraya-Watson estimator. Default \code{method="LL"}.}
  \item{tol}{ Tolerance parameter to avoid overflow when \code{nu} is larger than \code{tol}. Default is \code{tol=300}.}
  \item{from, to}{ Left and right-most points of the grid at which the density is to be estimated.}
  \item{len}{ Number of equally spaced points at which the density is to be estimated.}
}
\details{
  See Section 3 in Oliveira et al. (2013). See Di Marzio et al. (2009).
  The NAs will be automatically removed.
}
\value{
  Numeric vector of the same length of \code{t} with the values of the estimate at the evaluation points.
}
\references{ 
Di Marzio, M., Panzera A. and Taylor, C. C. (2009) Local polynomial regression for circular predictors. \emph{Statistics and Probability Letters}, 
\bold{79}, 2066--2075.

Oliveira, M., Crujeiras R.M. and Rodríguez--Casal, A. (2013) Nonparametric circular methods for exploring environmental data. 
\emph{Environmental and Ecological Statistics}, \bold{20}, 1--17..
}

\author{ María Oliveira, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\seealso{\code{\link{nu.LSCV.reg}}}

\examples{
data(speed.wind2)
dir <- rad(speed.wind2$Direction)
vel <- speed.wind2$Speed
nas <- which(is.na(vel))
dir <- dir[-nas]
vel <- vel[-nas]
t <- seq(0,2*pi,length=200)
estLL <- kern.reg.circ(dir, vel, t=t, nu=30)
estNW <- kern.reg.circ(dir, vel, t=t, nu=30, method="NW")

plot(dir, vel, xlab="direction", ylab="speed (m/s)", axes=FALSE)
lines(t, estLL, col=2)
lines(t, estNW, col=3)
axis(1, at=circular(seq(0,2*pi,by=pi/4)),
labels=c("N","NE","E","SE","S","SW","W","NW","N"))
legend("topleft", c("LL","NW"), lty=1, col=2:3)
axis(2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ circular regression }% __ONLY ONE__ keyword per line
