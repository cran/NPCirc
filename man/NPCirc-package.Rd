\encoding{utf8}
\name{NPCirc-package}
\alias{NPCirc-package}
\alias{NPCirc}
\docType{package}
\title{
Nonparametric circular methods.
}
\description{
Nonparametric smoothing methods for density and regression estimation involving circular data,
including methods for estimating circular densities and mean regression functions, described in Oliveira et al. (2014), testing proposals for circular mean regression described in Alonso-Pena et al. (2021), estimation of conditional characteristics (Alonso-Pena et al. 2022a,2022b) and estimation of conditional modes (Alonso-Pena and Crujeiras, 2022).
}
\details{
\tabular{ll}{
Package: \tab NPCirc\cr
Type: \tab Package\cr
Version: \tab 4.0.1\cr
Date: \tab 2021-07-21\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
This package incorporates the function \code{kern.den.circ} which computes the
circular kernel density estimator. For choosing the smoothing parameter different functions are available: \code{bw.rt}, \code{bw.CV},
\code{bw.pi}, and \code{bw.boot}. For mean regression involving circular variables, the package includes the functions: \code{kern.reg.circ.lin} for a circular
covariate and linear response; \code{kern.reg.circ.circ} for a circular covariate and a circular response; \code{kern.reg.lin.circ} for a linear covariate
and a circular response. The three functions compute Nadaraya-Watson and Local-Linear smoothers. The functions \code{bw.reg.circ.lin},
\code{bw.reg.circ.circ} and \code{bw.reg.circ.lin} implement cross--validation rules for selecting the smoothing parameter. Functions \code{circsizer.density} and \code{circsizer.re-} \code{gression} provide CircSiZer maps for kernel density estimation and regression estimation, respectively. Functions \code{noeffect.circ.lin}, \code{noeffect.circ.circ} and \code{noeffect.lin.circ} compute the test of no effect to assess the significance of the predictor variable. Additionally, functions \code{ancova.circ.lin}, \code{ancova.circ.circ} and \code{ancova.lin.circ} implement hypothesis testing tools to assess the equality and parallelism of regression curves across different groups of observations.

Function \code{circ.local.lik} implements the estimation of different functions of interest (transformations of the mean function) in contexts where the predictor variable is circular and the conditional distribution is a Gaussian, Bernoulli, Poisson or gamma. Function \code{bw.circ.local.lik} computes smoothing parameters for the estimation in the previously described cases, allowing for three different rules. Function \code{kern.dpreg.circ} implements the joint kernel estimation of the mean and dispersion functions in cases where the predictor is circular and the conditional distribution is a double Poisson, a particular case of the double exponential family. Smoothing parameters in this context can be computed with function \code{bw.joint.dpcirc}, using a two-step cross-validation algorithm.

Functions \code{modalreg.circ.circ}, \code{modalreg.lin.circ} and \code{modalreg.circ.lin} implement the estimation of the modal regression multifunction (the conditional local modes) in the three circular regression scenarios. Smoothing parameters for these contexts can by computed by modal cross-validation employin functions \code{bw.modalreg.circ.circ}, \code{bw.modalreg.lin.circ} and \code{bw.modalreg.circ.lin}.

Functions \code{dcircmix} and \code{rcircmix} compute the density function and generate random samples of a circular distribution or a mixture of circular
distributions, allowing for different components such as the circular uniform, von Mises, cardioid, wrapped Cauchy, wrapped normal and wrapped skew-normal.
Finally, some data sets are provided. Missing data are allowed. Registries with missing data are simply removed.

For a complete list of functions, use \code{library(help="NPCirc").} }
\author{
Maria Oliveira, Maria Alonso-Pena, Jose Ameijeiras-Alonso, Rosa M. Crujeiras, Alberto Rodriguez--Casal and Irene Gijbels.

Maintainer: Maria Alonso-Pena \email{mariaalonso.pena@usc.es}
}
\references{
Oliveira, M., Crujeiras, R.M. and Rodriguez--Casal, A. (2012) A plug--in rule for bandwidth selection in circular density. \emph{Computational Statistics
and Data Analysis}, \bold{56}, 3898--3908.

Oliveira, M., Crujeiras R.M. and Rodriguez--Casal, A. (2013) Nonparametric circular methods for exploring environmental data.
\emph{Environmental and Ecological Statistics}, \bold{20}, 1--17.

Oliveira, M., Crujeiras, R.M. and Rodriguez--Casal (2014) CircSiZer: an exploratory tool for circular data. \emph{Environmental and Ecological Statistics},
\bold{21}, 143--159.

Oliveira, M., Crujeiras R.M. and Rodriguez--Casal, A. (2014) NPCirc: an R package for nonparametric circular methods. \emph{Journal of Statistical Software},
\bold{61}(9), 1--26. \url{https://www.jstatsoft.org/v61/i09/}

Alonso-Pena, M., Ameijeiras-Alonso, J. and Crujeiras, R.M. (2021) Nonparametric tests for circular regression. \emph{Journal of Statistical Computation and Simulation}, \bold{91}(3), 477--500.

Alonso-Pena, M., Gijbels, I. and Crujeiras, R.M. (2022a). A general framework for circular local likelihood regression. \emph{Under review}.

Alonso-Pena, M., Gijbels, I. and Crujeiras, R.M. (2022b). Flexible joint modeling of mean and dispersion for the directional tuning of neuronal spike counts. \emph{Under review}.

Alonso-Pena, M. and Crujeiras, R. M. (2022). Analizing animal escape data with circular nonparametric multimodal regression. Annals of Applied Statistics. (To appear).
}



%\keyword{ circular density}
%\keyword{ circular regression}

\section{Acknowledgements}{ This work has been supported by Project MTM2008-03010  from the Spanish Ministry of Science; Project and MTM201676969-P from the AEI co-funded by the European Regional Development Fund (ERDF), the Competitive Reference
Groups 2017-2020 (ED431C 2017/38) from the Xunta de Galicia through the ERDF; and Innovation IAP network (Developing
crucial Statistical methods for Understanding major complex Dynamic Systems in natural, biomedical and social sciences (StUDyS)) from Belgian Science Policy. Work of Maria Alonso-Pena was supported by grant ED481A-2019/139 from Xunta de Galicia. Work of Jose Ameijeiras-Alonso was supported by the FWO research project G.0826.15N (Flemish Science Foundation); and GOA/12/014 project (Research Fund KU Leuven). The authors want to acknowledge Prof. Arthur Pewsey for facilitating data examples and for his comments.
}

