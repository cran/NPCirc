\encoding{latin1}
\name{bw.AA}
\alias{bw.AA}

\title{Direct and solve-the-equation plug-in rule}
\description{
Smoothing selectors for the circular kernel density (and its derivatives) estimator. This function implements the l-stage solve-the-equation and direct plug-in smoothing selector.
}
\usage{
bw.AA(x,deriv.order=0,method = c("ste","dpi"),nstage=2,kernel="vonmises",M=NULL,
      commonkappa=TRUE,Q1=NULL,Q2=NULL,lower=NULL,upper=NULL,tol=NULL,
      approximate=NULL)
}
\arguments{
  \item{x}{ Data from which the smoothing parameter is to be computed. The object is coerced to class \code{\link{circular}}.}
  \item{deriv.order}{ Derivative order. Default \code{deriv.order=0} (density estimation).}
  \item{method}{ Either \code{"ste"} (solve-the-equation) or \code{"dpi"} (direct plug-in). Can be abbreviated.}
  \item{nstage}{ Number of stages in the plug-in smoothing parameter. Default \code{nstage=2}.}
  \item{kernel}{ A character string giving the smoothing kernel to be used. This must be one of \code{vonmises} or \code{wrappednormal}. Default \code{kernel=vonmises}.}
  \item{M}{ Integer indicating the number of components in the von Mises mixture at stage 0. If \code{M} is a vector, AIC will be used, to select the number of components between the indicated values of \code{M}. Default, \code{M=1} if method \code{"ste"}, \code{M=1:5} if method \code{"dpi"}.}
  \item{commonkappa}{ Logical; if \code{TRUE}, at stage 0, all the components in the von Mises mixture have the same concentration. Default, \code{commonkappa=T}.}
  \item{Q1}{Vector of constants related to the kernel to derive the explicit expression of the optimal smoothing parameters of the density functionals. Its value is provided by default when using the \code{"vonmises"} or \code{"wrappednormal"} as \code{kernel}.}
  \item{Q2}{Constant related to the kernel to derive the explicit expression of the optimal smoothing parameters of the density derivative. Its value is provided by default when using the \code{"vonmises"} or \code{"wrappednormal"} as \code{kernel}.}
  \item{lower, upper}{ For method \code{"ste"}, the range over which the bandwidth h is searched. Default, \code{lower=10^(-3)} and \code{upper=pi^2/3} if \code{kernel="vonmises"}; \code{lower=10^(-10)} and \code{upper=0.5}, otherwise.}
  \item{tol}{ For method \code{"ste"}, the convergence tolerance for searching the smoothing parameter with \code{\link{uniroot}}. Default, \code{tol=.Machine$double.eps^0.25}.}
  \item{approximate}{ For method \code{"dpi"}, logical, if \code{TRUE}, the explicit expressions (relying on asymptotics) for the optimal smoothing parameters are employed. If \code{FALSE}, an optimization routine is employed, searching for the smoothing parameter minimizing the asymptotic mean squared error of the density derivatives and functionals. Default, \code{approximate=T}.}
}
\details{
By default, this function computes the solve-the-equation plug-in rule for circular kernel density estimation. If \code{method="dpi"}, this function computes the direct plug-in rule.

At stage 0, a mixture of von Mises is employed for computing the rule of thumb. The reason for employing a mixture model is that the von Mises estimates a uniform density when the true density model is k-fold rotational symmetric. Thus, in that case, the density functional estimator would be close to zero.

The number of components in the mixture is controlled with \code{M}. By default, a von Mises density (\code{M=1}) is employed in \code{method="ste"}. For \code{method="dpi"}, by default, the number of components in the mixture is selected using the Akaike Information Criterion, by comparing a mixture of 1 to 5 components. For simplicity, by default, the same concentration parameter is employed in all the components. This may be changed by setting \code{commonkappa=F}.


For \code{method="ste"}, the minimum number of stages is two (\code{nstage=2}). Otherwise, the solve-the-equation rule cannot be computed. A rule of thumb can be computed with \code{method="dpi"} and \code{nstage=0}.



}
\value{
Value of the smoothing parameter (mean resultant length). When the kernel is \code{vonmises} the bandwidth is equal to the concentration parameter.
}
\references{
Ameijeiras-Alonso, J. (2022). A reliable data-based smoothing parameter selection method for circular kernel estimation.
}
\author{
Jose Ameijeiras-Alonso.
}

\seealso{\code{\link{density.circular}}, \link{kern.den.circ}, \code{\link{bw.pi}}, \code{\link{bw.rt}}, \code{\link{bw.CV}}, \code{\link{bw.boot}} }
\examples{
set.seed(2022)
n <- 50
x <- rcircmix(n,model=6)
# Concentration parameter for density estimation
bw.AA(x) # Solve-the-equation concentration parameter
bw.AA(x, method="dpi") # Direct plug-in concentration parameter
# Concentration parameter for the density derivative estimate
bw.AA(x, method="ste") # Solve-the-equation concentration parameter
bw.AA(x, deriv.order=1, method="dpi") # Direct plug-in concentration parameter
}

\keyword{ circular density }
