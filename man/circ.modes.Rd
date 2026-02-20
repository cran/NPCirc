\name{circ.modes}
\alias{circ.modes}
\title{
Estimation of the modal set}
\description{
Computes the estimator of the modes of a circular variable based on a kernel density estimator. The estimation is performed numerically through the circular mean shift algorithm.
}
\usage{
circ.modes(x,bw.modes=NULL,bw.dens=NULL,plot=TRUE,col.modes=2,col.dens=1,lwd.modes=2,
            lwd.dens=2,tol=0.0000000000001, labels=NULL, control.circular = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data from which the modes are to be estimated. The object is coerced to class \code{\link[circular]{circular}}.}
  \item{bw.modes}{Smoothing parameter for the von Mises kernel for the estimation of the modes. If null, \code{bw.AA} with the option \code{deriv=1} is used.}
  \item{bw.dens}{Smoothing parameter for the von Mises kernel for the estimation of the density. If null, \code{bw.AA} with the option \code{deriv=0} is used. Only used if \code{plot = TRUE}.}
  \item{plot}{Logical; if \code{TRUE}, a plot of the estimated density and detected modes is produced.}
  \item{col.modes}{Color for the mode lines. Can be specified as a character string (color name or hexadecimal code) or as a numeric index. Only used if \code{plot = TRUE}.}
  \item{col.dens}{Color for the density curve. Can be specified as a character string (color name or hexadecimal code) or as a numeric index. Only used if \code{plot = TRUE}.}
  \item{lwd.modes}{Line width for the mode lines (see graphical parameter \code{lwd}). Must be a positive numeric value. Only used if \code{plot = TRUE}.}
  \item{lwd.dens}{Line width for the density curve (see graphical parameter \code{lwd}). Must be a positive numeric value. Only used if \code{plot = TRUE}.}
  \item{tol}{Numeric tolerance used in the mode estimation algorithm to determine convergence. Smaller values yield more precise but potentially slower computations.}
  \item{labels}{Logical; if \code{TRUE}, angular reference labels are added to the plot.}
  \item{control.circular}{A list of control parameters passed to \code{circular::circular} to define properties such as units, zero direction, and rotation.}
}

\details{
The function performs nonparametric density estimation for circular or axial data using kernel methods and identifies modes of the estimated density. Mode detection is carried out through an iterative procedure that searches for local maxima of the smoothed density function.

The number and location of the detected modes depend on the smoothing parameter used in the density estimation. In particular, larger bandwidth values produce smoother densities with fewer modes, whereas smaller bandwidth values may reveal additional local features. If \code{bw=NULL}, the bandwidth is selected as a plug-in of the optimal bandwidth for the estimation of the first derivative. 

When \code{plot = TRUE}, the estimated density curve and the detected modes are displayed. Graphical parameters such as colors and line widths can be customized through the corresponding arguments.

}


\value{
A list with the following components:

\item{data}{The input data.}

\item{modes}{Numeric vector containing the estimated mode locations (in the same units as the input data).}

\item{bw.modes}{Smoothing parameter used for mode estimation.}

\item{bw.dens}{Smoothing parameter used for density estimation.}
}

\references{
Alonso-Pena, M., Aerden, D.G.A.M. and Angulo, J.M. (2026). Quantifying tectonic relevance of porphyroblast inclusion trails: axial mode estimation and clustering. 
}
\author{
Maria Alonso-Pena
}




\examples{
# axial data in degrees, geographics template
set.seed(1)
x<-c(rvonmises(100,pi/8,10),rvonmises(100,pi,5))
x<-conversion.circular(x/2,units="degrees",template="geographics",modulo="pi")

plot(x,stack=TRUE)
points(x+180,stack=TRUE,col="grey")

circ.modes(x,col.modes=3,col.dens="blue",lwd.dens=3,lwd.modes=3)

}

\keyword{ mode estimation }
\keyword{ circular mean shift  }

