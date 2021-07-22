\name{noeffect.circ.lin}
\alias{noeffect.circ.lin}
\alias{noeffect.lin.circ}
\alias{noeffect.circ.circ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{No-effect test for regression with circular data}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function \code{noeffect.circ.lin} computes the no-effect test for a circular predictor variable and a real-valued response variable as described in Alonso-Pena et al. (2021). It uses the nonparametric Nadaraya-Watson estimator or the Local-Linear estimator for circular-linear data described in Di Marzio et al. (2009) and Oliveira et al. (2013).

Function \code{noeffect.lin.circ} computes the no-effect test for a real-valued predictor variable and a circular response variable as described in Alonso-Pena et al. (2021). It uses the nonparametric Nadaraya-Watson estimator or the Local-Linear estimator for linear-circular data described in Di Marzio et al. (2012).

Function \code{noeffect.circ.circ} computes the no-effect test for a circular predictor variable and a circular response variable as described in Alonso-Pena et al. (2021). It uses the nonparametric Nadaraya-Watson estimator or the Local-Linear estimator for circular-circular data described in Di Marzio et al. (2012).
}
\usage{
noeffect.circ.lin(x, y, bw, method = "LL", calib = "chisq", n_boot = 500)
noeffect.lin.circ(x, y, bw, method = "LL", n_boot = 500)
noeffect.circ.circ(x, y, bw, method = "LL", n_boot = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of data for the independent variable. The object is coerced to              class \code{\link{circular}} when using functions \code{noeffect.circ.lin} and \code{noeffect.circ.circ}.}
  \item{y}{Vector of data for the dependent variable. This must be same length as             \code{x}. The object is coerced to class \code{\link{circular}} when using functions \code{noeffect.lin.circ} and \code{noeffect.circ.circ}.}
  \item{bw}{Smoothing parameter to be used. If not provided, functions \code{noeffect.circ.lin} and  \code{noeffect.circ.circ} select \eqn{4cv} and function  \code{noeffect.lin.circ} selects \eqn{cv/4},  where \eqn{cv} is the parameter selected by cross-validation.}
  \item{method}{Character string giving the estimator to be used. This must be one of \code{"LL"} for Local-Linear estimator or \code{"NW"} for Nadaraya-Watson estimator. Default \code{method="LL"}.}
  \item{calib}{Character string giving the calibration method to be used in \code{noeffect.circ.lin} function. This must be one of \code{"chisq"} for the chi-squared approximation or \code{"boot"} for the bootstrap calibration.}
  \item{n_boot}{Number of bootstrap resamples. Default is \code{n_boot=500}. In function \code{noeffect.circ.lin}, only if \code{calib="boot"}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
See Alonso-Pena et al. (2021). The NAs will be automatically removed.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list with class \code{"htest"} containing the following components:
  \item{statistic}{observed value of the statistic.}
  \item{bw}{Smoothing parameter used.}
  \item{p.value}{p-value for the test.}
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{alternative}{a character string describing the alternative hypothesis.}
}
\references{
%% ~put references to the literature/web site here ~
Alonso-Pena, M., Ameijeiras-Alonso, J. and Crujeiras, R.M. (2021) Nonparametric tests for circular regression. \emph{Journal of Statistical Computation and Simulation}, \strong{91}(3), 477--500.

Di Marzio, M., Panzera A. and Taylor, C. C. (2009) Local polynomial regression for circular predictors. \emph{Statistics and Probability Letters}, \strong{79}, 2066--2075.

Di Marzio, M., Panzera A. and Taylor, C. C. (2012) Non--parametric regression for circular responses. \emph{Scandinavian Journal of Statistics}, \strong{40}, 228--255.

Oliveira, M., Crujeiras R.M. and Rodriguez-Casal, A. (2013) Nonparametric circular methods for exploring environmental data. \emph{Environmental and Ecological Statistics}, \strong{20}, 1--17.


}
\author{
%%  ~~who you are~~
María Alonso-Pena, Jose Ameijeiras-Alonso and Rosa M. Crujeiras
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{kern.reg.circ.lin}}, \code{\link{kern.reg.lin.circ}}, \code{\link{kern.reg.circ.circ}}
}
\examples{
\donttest{
# No-effect circ-lin
set.seed(2025)
x <- rcircularuniform(200)
y <- 2*sin(as.numeric(x)) + rnorm(200, sd=2)
noeffect.circ.lin(x, y)

# No-effect lin-circ
set.seed(2025)
x <- runif(200)
y<- pi/8 + rvonmises(200, mu = 0, kappa = 0.75)
noeffect.lin.circ(x, y)


# No-effect circ-circ
set.seed(2025)
x <- rcircularuniform(200)
y <- atan2(sin(2*x),cos(2*x)) + rvonmises(200, mu = 0, kappa = 2)
noeffect.circ.circ(x, y)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~htest }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~smooth }% __ONLY ONE__ keyword per line
