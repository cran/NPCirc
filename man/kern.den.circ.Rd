\encoding{latin1}
\name{kern.den.circ}
\alias{kern.den.circ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Nonparametric circular kernel density estimation }
\description{
   This function computes circular kernel estimates with the given bandwidth, taking the von Mises distribution as circular kernel.
}
\usage{
kern.den.circ(x, t=NULL, nu, from=0, to=2*pi, len=250)
}
\arguments{
 
  \item{x}{ Sample of angles in radians (between \eqn{0} and \eqn{2\pi}) from which the estimate is to be computed.}
  \item{t}{ Vector of angles in radians where to evaluate the estimator. If \code{NULL} equally spaced points are used according to the parameters 
            \code{from}, \code{to} and \code{len}.}
  \item{nu}{ Smoothing parameter to be used. The value of \code{nu} can be chosen by using the functions \code{\link{nu.rt}}, \code{\link{nu.CV}}, 
             \code{\link{nu.pi}} and \code{\link{nu.boot}}.}
  \item{from, to}{ Left and right-most points of the grid at which the density is to be estimated.}
  \item{len}{ Number of equally spaced points at which the density is to be estimated.}
}
\details{
  The NAs will be automatically removed.
}
\value{
  Numeric vector of the same length of \code{t} with the kernel density estimated values at \code{t}.
}
\references{ 
Oliveira, M., Crujeiras, R.M. and Rodríguez--Casal, A. (2012) A plug--in rule for bandwidth selection in circular density. \emph{Computational Statistics 
and Data Analysis}, \bold{56}, 3898--3908.

Taylor, C.C. (2008) Automatic bandwidth selection for circular density estimation. \emph{Computational Statistics and Data Analysis}, \bold{52}, 3493--3500.
}


\author{ María Oliveira, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\seealso{ \code{\link{nu.rt}}, \code{\link{nu.CV}}, \code{\link{nu.pi}}, \code{\link{nu.boot}} }

\examples{
##  Estimating the density function of a sample of circular data
set.seed(2012)
n <- 100
x <- rcircmix(n, model=14)
t <- seq(0,2*pi,length=500)
est <- kern.den.circ(x, t, nu=50)
plot(t, dcircmix(t,model=14), ylim=c(0,0.4), type="l", lwd=2, main="Linear representation")
lines(t, est, col=2)
plot(circular(x), shrink=1.3, main="Circular representation")
lines(circular(t), dcircmix(t,model=14), lwd=2)
lines(circular(t), est, col=2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ circular density }% __ONLY ONE__ keyword per line
