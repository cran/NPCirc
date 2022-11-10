\name{ancova.circ.lin}
\alias{ancova.circ.lin}
\alias{ancova.lin.circ}
\alias{ancova.circ.circ}
\title{Nonparametric analysis of covariance tests for circular regression}
\description{
Function \code{ancova.circ.lin} computes nonparametric ANCOVA tests to compare regression curves with a circular predictor variable and a real-valued response variable. The null hypothesis may be either equality or parallelism of the regression curves, as described in Alonso-Pena et al. (2021). It uses the nonparametric Nadaraya-Watson estimator or the Local-Linear estimator for circular-linear data described in Di Marzio et al. (2009) and Oliveira et al. (2013).

Function \code{ancova.lin.circ} computes nonparametric ANCOVA tests to compare regression curves with a real-valued predictor variable and a circular response variable. The null hypothesis may be either equality or parallelism of the regression curves, as described in Alonso-Pena et al. (2021). It uses the nonparametric Nadaraya-Watson estimator or the Local-Linear estimator for linear-circular data described in Di Marzio et al. (2012).

Function \code{ancova.circ.circ} computes nonparametric ANCOVA tests to compare regression curves with a circular predictor variable and a circular response variable. The null hypothesis may be either equality or parallelism of the regression curves, as described in Alonso-Pena et al. (2021). It uses the nonparametric Nadaraya-Watson estimator or the Local-Linear estimator for circular-circular data described in Di Marzio et al. (2012).
}
\usage{
ancova.circ.lin(x, y, g, bw, bw1, test = "eq", method = "LL",calib = "chisq", 
   n_boot = 500)
ancova.lin.circ(x, y, g, bw, bw1, test = "eq", method = "LL", n_boot = 500)
ancova.circ.circ(x, y, g, bw, bw1, test = "eq", method = "LL", n_boot = 500)
}
\arguments{
\item{x}{Vector of data for the independent variable. The object is coerced to class \code{\link{circular}} when using functions \code{ancova.circ.lin} and \code{ancova.circ.circ}.}
\item{y}{Vector of data for the dependent variable. This must be same length as \code{x}. The object is coerced to class \code{\link{circular}} when using functions \code{ancova.lin.circ} and \code{ancova.circ.circ}.}
\item{g}{Vector of group indicators.}
\item{bw}{Smoothing parameter to be used. If not provided it selects the parameter obtained by cross-validation.}
\item{bw1}{Preliminary smoothing parameter for the parallelism test.}
\item{test}{Character string giving the type of test to be performed. Must be one of \code{"eq"} for the test of equality or \code{"paral"} for the test of parallelism. Default is \code{test="eq"}.}
\item{method}{Character string giving the estimator to be used. This must be one of \code{"LL"} for Local-Linear estimator or \code{"NW"} for Nadaraya-Watson estimator. Default \code{method="LL"}.}
\item{calib}{Character string giving the calibration method to be used in \code{ancova.circ.lin} function. This must be one of \code{"chisq"} for the chi-squared approximation or \code{"boot"} for the bootstrap calibration.}
\item{n_boot}{Number of bootstrap resamples. Default is \code{n_boot=500}. In function \code{ancova.circ.lin}, only if \code{calib="boot"}.}
}
\details{
See Alonso-Pena et al. (2021). The NAs will be automatically removed.
}
\value{
A list with class \code{"htest"} containing the following components:
\item{statistic}{observed value of the statistic.}
\item{bw}{Smoothing parameter used.}
\item{p.value}{p-value for the test.}
\item{data.name}{a character string giving the name(s) of the data.}
\item{alternative}{a character string describing the alternative hypothesis.}
}
\references{
Alonso-Pena, M., Ameijeiras-Alonso, J. and Crujeiras, R.M. (2021) Nonparametric tests for circular regression. \emph{Journal of Statistical Computation and Simulation}, \bold{91}(3), 477--500.

Di Marzio, M., Panzera A. and Taylor, C. C. (2009) Local polynomial regression for circular predictors. \emph{Statistics and Probability Letters}, \bold{79}, 2066--2075.

Di Marzio, M., Panzera A. and Taylor, C. C. (2012) Non--parametric regression for circular responses. \emph{Scandinavian Journal of Statistics}, \bold{40}, 228--255.

Oliveira, M., Crujeiras R.M. and Rodriguez-Casal, A. (2013) Nonparametric circular methods for exploring environmental data. \emph{Environmental and Ecological Statistics}, \bold{20}, 1--17.
}
\author{
Maria Alonso-Pena, Jose Ameijeiras-Alonso and Rosa M. Crujeiras
}
\seealso{
\code{\link{kern.reg.circ.lin}}, \code{\link{kern.reg.lin.circ}}, \code{\link{kern.reg.circ.circ}}
}
\examples{
\donttest{
# ANCOVA circ-lin
set.seed(2025)
x1 <- rcircularuniform(100)
x2 <- rcircularuniform(100)
x <- c(x1, x2)
y1 <- 2*sin(as.numeric(x1)) + rnorm(100, sd=2)
y2 <- 4 + 2*sin(as.numeric(x2)) + rnorm(100, sd=2)
y <- c(y1, y2)
g <- c(rep(0,100), rep(1,100))
ancova.circ.lin(x, y, g,  test = "eq")
ancova.circ.lin(x, y, g,  test = "paral")


# ANCOVA lin-circ
set.seed(2025)
x1 <- runif(100)
x2 <- runif(100)
y1 <- 3*pi*x1^2 + rvonmises(100, mu = 0, kappa = 6)
y2 <- 2*pi/8 + 3*pi*x2^2 + rvonmises(100, mu = 0, kappa = 6)
x <- c(x1, x2)
y <- c(y1, y2)
g<-c(rep(0, 100), rep(1, 100))
ancova.lin.circ(x, y, g,  test = "eq")
ancova.lin.circ(x, y, g,  test = "paral")


# ANCOVA circ-circ
set.seed(2025)
x1 <- rcircularuniform(100)
x2 <- rcircularuniform(100)
y1 <- 2*sin(2*x1) + rvonmises(100, mu = 0, kappa = 8 )
y2 <- pi/8 + 2*sin(2*x2) + rvonmises(100, mu = 0, kappa = 8 )
x <- c(x1, x2)
y <- c(y1, y2)
g<-c(rep(0, 100), rep(1, 100))
ancova.circ.circ(x, y, g,  test = "eq")
ancova.circ.circ(x, y, g,  test = "paral")
}
}
\keyword{ ~htest }
\keyword{ ~smooth }