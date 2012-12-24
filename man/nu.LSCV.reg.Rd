\encoding{latin1}
\name{nu.LSCV.reg}
\alias{nu.LSCV.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Least squares cross-validation for circular-linear regression estimation }
\description{
  This function provides the least squares cross-validation (LSCV) smoothing parameter for the Nadaraya-Watson and Local-Linear estimators when the covariate
  is circular and the response variable is linear.
}
\usage{
nu.LSCV.reg(x, y, method="LL", lower=0, upper=100, tol=0.1)
}
\arguments{
  \item{x}{ Vector of data in radians (between \eqn{0} and \eqn{2\pi}) for the circular covariate.}
  \item{y}{ Vector of data for the dependent linear variable. This must be same length as \code{x}.}
  \item{method}{ Character string giving the estimator to be used. This must be one of \code{"LL"} or \code{"NW"}. Default \code{method="LL"}.}
  \item{lower, upper}{ \code{lower} and \code{upper} boundary of the interval to be used in the search for the value of the smoothing parameter. Default
                       \code{lower=0} and \code{upper=100}.}
  \item{tol}{ Convergence tolerance for \code{\link{optimize}}. Default \code{tol=0.1}.}
}
\details{
  The LSCV smoothing parameter is obtained as the value of \eqn{\nu} that minimizes expression (13) in Oliveira et al. (2013). 
  The NAs will be automatically removed.
}
\value{
  Value of the smoothing parameter.
}
\references{ 
Oliveira, M., Crujeiras R.M. and Rodríguez--Casal, A. (2013) Nonparametric circular methods for exploring environmental data. 
\emph{Environmental and Ecological Statistics}, \bold{20}, 1--17.
}

\author{ María Oliveira, Rosa M. Crujeiras and Alberto Rodríguez--Casal}

\seealso{\code{\link{kern.reg.circ}} }

\examples{
set.seed(2012)
n <- 100
x <- seq(0,2*pi,length=n)
y <- sin(x)+0.2*rnorm(n)
nu.LSCV.reg(x, y, method="LL", lower=1, upper=20)
nu.LSCV.reg(x, y, method="NW", lower=1, upper=20)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ circular regression }% __ONLY ONE__ keyword per line
