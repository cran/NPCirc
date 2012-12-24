\encoding{latin1}
\name{NPCirc-package}
\alias{NPCirc-package}
\alias{NPCirc}
\docType{package}
\title{
Nonparametric circular methods.
}
\description{
This package implements nonparametric kernel methods for density and regression estimation for circular data.
}
\details{
\tabular{ll}{
Package: \tab NPCirc\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2012-12-24\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
This package incorporates the function \code{kern.den.circ} which computes the 
kernel circular density estimator. For choosing the smoothing parameter different functions are available: \code{nu.rt}, \code{nu.CV},
\code{nu.pi}, and \code{nu.boot}. For circular-linear regression (circular covariate and linear response), Nadaraya-Watson and 
Local-Linear smoothers, are also available in function \code{kern.reg.circ}. The function \code{nu.LSCV.reg} computes the least squares
cross-validation bandwidth for those estimators. 
Functions \code{circsizer.density} and \code{circsizer.regression} provides CircSiZer maps for kernel density estimation and regression estimation, respectively. 
Functions \code{dcircmix} and \code{rcircmix} compute the density function and generate random samples of a circular distribution or a mixture of circular 
distributions, allowing for different components such as the circular uniform, von Mises, cardioid, wrapped Cauchy, wrapped normal and wrapped skew-normal. 
Finally, some data sets are provided. Missing data are allowed. Registries with missing data are simply removed. 

For a complete list of functions, use \code{library(help="NPCirc").} }
\author{
María Oliveira, Rosa M. Crujeiras and Alberto Rodríguez--Casal

Maintainer: María Oliveira \email{maria.oliveira@usc.es}
}
\references{
Oliveira, M., Crujeiras, R.M. and Rodríguez--Casal, A. (2012) A plug--in rule for bandwidth selection in circular density. \emph{Computational Statistics 
and Data Analysis}, \bold{56}, 3898--3908.

Oliveira, M., Crujeiras R.M. and Rodríguez--Casal, A. (2013) Nonparametric circular methods for exploring environmental data. 
\emph{Environmental and Ecological Statistics}, \bold{20}, 1--17.

Oliveira, M., Crujeiras, R.M. and Rodríguez--Casal (2013) CircSiZer: an exploratory tool for circular data. http://arxiv.org/abs/1210.2555v2.
}


%\keyword{ circular density}
%\keyword{ circular regression}

\section{Acknowledgements}{ This work has been supported by Project MTM2008-03010 from the Spanish Ministry of Science and Innovation IAP network (Developing 
crucial Statistical methods for Understanding major complex Dynamic Systems in natural, biomedical and social sciences (StUDyS)) from Belgian Science Policy.
The authors want to acknowledge Prof. Arthur Pewsey for facilitating data examples and for his comments.
}

