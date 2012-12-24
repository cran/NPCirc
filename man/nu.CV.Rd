\encoding{latin1}
\name{nu.CV}
\alias{nu.CV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cross-validation for density estimation}
\description{
  This function gives the least squares cross-validation smoothing parameter or the likelihood cross-validation smoothing parameter for density estimation.
}
\usage{
nu.CV(x, method="LCV", lower=0, upper=100, tol=0.1, np=500)
}
\arguments{
  \item{x}{ Sample of angles in radians (between \eqn{0} and \eqn{2\pi}) from which the smoothing parameter is to be computed.}
  \item{method}{ Character string giving the cross-validation rule to be used. This must be one of \code{"LCV"} or \code{"LSCV"}. Default \code{method="LCV"}.}
  \item{lower, upper}{ \code{lower} and \code{upper} boundary of the interval to be used in the search for the value of the smoothing parameter. Default
                       \code{lower=0} and \code{upper=100}.}
  \item{tol}{ Convergence tolerance for \code{\link{optimize}}. Default \code{tol=0.1}.}
  \item{np}{  Number of points where to evaluate the estimator for numerical integration when \code{method="LSCV"}. Default \code{np=500}.}
}
\details{
  The LCV smoothing parameter is obtained as the value of \eqn{\nu} that maximizes the logarithm of the likelihood cross-validation function (8) in Oliveira 
  et al. (2013). 
  The LSCV smoothing parameter is obtained as the value of \eqn{\nu} that minimizes expression (7) in Oliveira et al. (2013). 
  See also Hall et al. (1987) and Oliveira et al. (2012).
  The NAs will be automatically removed.
}
\value{
  Value of the smoothing parameter.
}
\references{ 
Hall, P., Watson, G.S. and Cabrera, J. (1987) Kernel density estimation with spherical data, \emph{Biometrika}, \bold{74}, 751--762.

Oliveira, M., Crujeiras, R.M. and Rodríguez--Casal, A. (2012) A plug--in rule for bandwidth selection in circular density. \emph{Computational Statistics 
and Data Analysis}, \bold{56}, 3898--3908.

Oliveira, M., Crujeiras R.M. and Rodríguez--Casal, A. (2013) Nonparametric circular methods for exploring environmental data. 
\emph{Environmental and Ecological Statistics}, \bold{20}, 1--17.
}

\seealso{\code{\link{kern.den.circ}}, \code{\link{nu.rt}}, \code{\link{nu.pi}}, \code{\link{nu.boot}} }

\author{ María Oliveira, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\examples{
set.seed(2012)
n <- 100
x <- rcircmix(n, model=11)
nu.CV(x, method="LCV", lower=0, upper=20)
nu.CV(x, method="LSCV", lower=0, upper=20)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ circular density }% __ONLY ONE__ keyword per line
