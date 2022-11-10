\name{kern.den.circ}
\alias{kern.den.circ}

\title{
Kernel density derivative estimate for circular data
}
\description{
This function computes the kernel density derivative estimate with the given kernel and bandwidth for circular data.
}
\usage{
kern.den.circ(x,z=NULL,bw="AA",deriv.order=0,kernel="vonmises",na.rm = FALSE,
              from = circular(0), to = circular(2 * pi),n = 512,control.circular=list())
}

\arguments{
  \item{x}{ Data from which the estimate is to be computed. The object is coerced to class \code{\link{circular}}.}
  \item{z}{ Points where the density derivative is estimated. If \code{NULL}
    equally spaced points are used according to the parameters
    \code{from}, \code{to} and \code{n}.}
  \item{bw}{ Smoothing parameter to be used. \code{bw} can also be a character string giving a rule to choose the bandwidth. See \code{\link{bw.AA}}, \code{\link{bw.pi}}, \code{\link{bw.rt}}, \code{\link{bw.CV}}, and \code{\link{bw.boot}}. The default, \code{'AA'}, is the 2-stage solve-the-equation plug-in smoothing selector.}
  \item{deriv.order}{ Derivative order. Default \code{deriv.order=0} (density estimation).}
  \item{kernel}{a character string giving the smoothing kernel to be
          used. This must be one of \code{"vonmises"} or
      \code{"wrappednormal"}.}
  \item{na.rm}{logical; if \code{TRUE}, missing values are removed from
    \code{x}. If \code{FALSE} any missing values cause an error.}
  \item{from, to}{the left and right-most
    points of the grid at which the density is to be estimated. The objects are coerced to class \code{circular}.}
  \item{n}{the number of equally spaced points at which the density is
          to be estimated.}
  \item{control.circular}{the attribute of the resulting objects (\code{x} component).}
}

\value{
  An object with class \code{\link{density.circular}} whose
  underlying structure is a list containing the following components.
  \item{data}{original dataset.}
  \item{x}{the \code{n} coordinates of the points where the density is
    estimated. It is a circular objects with coordinate system setting using \code{control.circular}.}
  \item{y}{the estimated density values.}
  \item{bw}{the smoothing parameter used.}
  \item{N}{the sample size after elimination of missing values.}
  \item{call}{the call which produced the result.}
  \item{data.name}{the deparsed name of the \code{x} argument.}
  \item{has.na}{logical, for compatibility (always FALSE).}
}
\references{
Ameijeiras-Alonso, J. (2022) A reliable data-based smoothing parameter selection method for circular kernel estimation.

Di Marzio, M., Panzera, A., & Taylor, C. C. (2011). Kernel density estimation on the torus. \emph{Journal of Statistical Planning and Inference}, \bold{141}(6), 2156--2173.

Oliveira, M., Crujeiras R.M. and Rodr?guez-Casal, A. (2014) NPCirc: an R package for nonparametric circular methods. \emph{Journal of Statistical Software},
\bold{61}(9), 1--26. \url{https://www.jstatsoft.org/v61/i09/}
}
\author{
Jose Ameijeiras-Alonso.
}

\seealso{\code{\link{bw.AA}}, \code{\link{plot.density.circular}}, \code{\link{lines.density.circular}}, \code{\link{bw.pi}}, \code{\link{bw.rt}}, \code{\link{bw.CV}}, \code{\link{bw.boot}} }
\examples{
set.seed(2022)
n <- 50
x <- rcircmix(n, model=13)
# Using the smoothing parameter by default,
# i.e., 2-stage solve-the-equation plug-in rule
est1 <- kern.den.circ(x,deriv.order=1)
# Selecting the smoothing parameter: 2-stage direct plug-in rule
est2 <- kern.den.circ(x, bw="dpi", deriv.order=1)
# Circular plot
plot(est1, plot.type="circle", points.plot=TRUE, shrink=1.4,
      main="Circular plot",ylab="Density derivative circular")
lines(est2, plot.type="circle", shrink=1.4 ,col=2)
# Linear plot
plot(est1, plot.type="line", main="Linear plot",ylab="Density derivative circular")
lines(est2, plot.type="line", col=2)
}

\keyword{ circular density }
