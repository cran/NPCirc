\name{modal.clustering}
\alias{modal.clustering}
\title{
Unsupervised clustering for circular or axial data}
\description{
Performs nonparametric modal clustering for circular or axial data using kernel density estimation. Observations are assigned to clusters according to the modes of the estimated density, without requiring the number of clusters to be specified in advance.
}
\usage{
modal.clustering(x,bw=NULL,tol=0.0000000000001, plot=TRUE,shrink = 1, pch=16, 
      col.points="grey",stack=TRUE, labels=NULL, control.circular = list())
}
\arguments{
  \item{x}{Data from which the modes are to be estimated. The object is coerced to class \code{\link[circular]{circular}}.}
  \item{bw}{Smoothing parameter for the von Mises kernel for the estimation of the modes. If null, \code{bw.AA} with the option \code{deriv=1} is used.}
  \item{tol}{Numeric tolerance used in the mode estimation algorithm to determine convergence. Smaller values yield more precise but potentially slower computations.}
  \item{plot}{Logical; if \code{TRUE}, a plot of the estimated density and detected modes is produced.}
  \item{shrink}{Positive numeric value controlling the radial scaling of stacked points in the plot. Only used if \code{plot = TRUE}.}
  \item{pch}{Plotting character for the data points (see graphical parameter \code{pch}). Only used if \code{plot = TRUE}.}
  \item{col.points}{Color of the data points. Can be specified as a character string or numeric index. Only used if \code{plot = TRUE}.}
  \item{stack}{Logical; if \code{TRUE}, observations are stacked radially in the plot to reduce overlap. Only used if \code{plot = TRUE}.}
  \item{labels}{Logical; if \code{TRUE}, angular reference labels are added to the plot.}
  \item{control.circular}{A list of control parameters passed to \code{circular::circular} to define properties such as units, zero direction, and rotation.}
}

\details{
The function performs nonparametric mode-based clustering for circular or axial data using kernel density estimation. The circular mean-shift algorithm is initialized with all sample points as starting values. Points that converge to the same mode are assigned to the same cluster, providing a fully data-driven clustering solution without requiring the number of clusters to be specified in advance.

The bandwidth parameter (\code{bw}) controls the smoothness of the estimated density and strongly influences the number of detected modes. If \code{bw=NULL}, the bandwidth is selected as a plug-in of the optimal bandwidth for the estimation of the first derivative. The \code{tol} argument defines the numerical tolerance for convergence in the mode-seeking algorithm.

When \code{plot = TRUE}, the estimated density and detected modes are displayed in a circular plot. Data points can be optionally stacked radially (\code{stack}) to reduce overlap, and their appearance can be customized through \code{pch}, \code{col.points}, and \code{shrink}. Angular reference labels can be added using \code{labels}. Only the plotting-related arguments are used when \code{plot = TRUE}.

The \code{control.circular} list can be used to specify properties of the circular or axial data, such as units, zero direction, and rotation, ensuring that calculations respect the geometry of the sample.
}


\value{
A list with the following components:

\item{modes}{Numeric vector containing the estimated mode locations (in the same units as the input data).}

\item{bw}{Smoothing parameter used for mode estimation.}

\item{clust}{An integer vector assigning each observation in \code{x} to a cluster. Observations assigned the same number belong to the same cluster.}
}

\references{
Alonso-Pena, M., Aerden, D.G.A.M. and Angulo, J.M. (2026). Quantifying tectonic relevance of porphyroblast inclusion trails: axial mode estimation and clustering. 
}
\author{
Maria Alonso-Pena
}




\examples{
\donttest{
# axial data in degrees, geographics template
set.seed(1)
x<-c(rvonmises(100,pi/8,10),rvonmises(100,pi,5))
x<-conversion.circular(x/2,units="degrees",template="geographics",modulo="pi")

plot(x,stack=TRUE)
points(x+180,stack=TRUE,col="grey")

modal.clustering(x,col.points="grey80")

}
}

\keyword{ unsupervised clustering }
\keyword{ circular mean shift  }

