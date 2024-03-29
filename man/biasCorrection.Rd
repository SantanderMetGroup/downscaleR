% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biasCorrection.R
\name{biasCorrection}
\alias{biasCorrection}
\title{Bias correction methods}
\usage{
biasCorrection(
  y,
  x,
  newdata = NULL,
  precipitation = FALSE,
  method = c("delta", "scaling", "eqm", "pqm", "gpqm", "loci", "dqm", "qdm", "isimip3"),
  cross.val = c("none", "loo", "kfold"),
  folds = NULL,
  consecutive = TRUE,
  window = NULL,
  scaling.type = c("additive", "multiplicative"),
  fitdistr.args = list(densfun = "normal"),
  wet.threshold = 1,
  n.quantiles = NULL,
  extrapolation = c("none", "constant"),
  theta = c(0.95, 0.05),
  detrend = TRUE,
  isimip3.args = NULL,
  join.members = FALSE,
  return.raw = FALSE,
  interpGrid.args = list(),
  parallel = FALSE,
  max.ncores = 16,
  ncores = NULL
)
}
\arguments{
\item{y}{A grid or station data containing the observed climate data for the training period}

\item{x}{A grid containing the simulated climate by the model for the training period. This can be either
the same variable as \code{obs}, in the case of model calibration (bias correction and related techniques) or MOS (model
 output statistics) approaches, or a set of predictors in the case of \emph{perfect prog} downscaling approaches (possibly
 after principal component analysis via the \code{\link{prinComp}} output).}

\item{newdata}{A grid containing the simulated climate for the variables used in \code{x}, but considering the test period.
If \code{x} is a multimember grid (the member dimension is > 1), the function assumes that the order of 
the members is the same in \code{newdata}.}

\item{precipitation}{Logical for precipitation data (default to FALSE). If TRUE adjusts precipitation 
frequency in 'x' (prediction) to the observed frequency in 'y' (see Details). To adjust the frequency, 
parameter \code{wet.threshold} is used (see below).}

\item{method}{method applied. Current accepted values are \code{"eqm"}, \code{"delta"},
\code{"scaling"}, \code{"pqm"} and \code{"gpqm"} \code{"variance"},\code{"loci"}, \code{"ptr"}, 
\code{"dqm"}, \code{"qdm"}, \code{"isimip3"}. See details.}

\item{cross.val}{Logical (default to FALSE). Should cross-validation be performed? methods available are 
leave-one-out ("loo") and k-fold ("kfold") on an annual basis. The default option ("none") does not 
perform cross-validation.}

\item{folds}{Only requiered if \code{cross.val = "kfold"}. Integer indicating the number of folds (see 
argument \code{consecutive}) or a list of vectors, each containing the years to be grouped in the corresponding fold.}

\item{consecutive}{Default is TRUE. Create folds containing consecutive years? Only used if cross.val = "kfold" and
folds is an integer. If FALSE, each years will be sampled randomly to create the folds.}

\item{window}{vector of length = 2 (or 1) specifying the time window width used to calibrate and the 
target days (days that are being corrected). The window is centered on the target day/s 
(window width >= target days). Default to \code{NULL}, which considers the whole period.}

\item{scaling.type}{Character indicating the type of the scaling method. Options are \code{"additive"} (default)
or \code{"multiplicative"} (see details). This argument is ignored if \code{"scaling"} is not 
selected as the bias correction method.}

\item{fitdistr.args}{Further arguments passed to function \code{\link[MASS]{fitdistr}} 
(\code{densfun}, \code{start}, \code{...}). Only used when applying the "pqm" method 
(parametric quantile mapping). Please, read the \code{\link[MASS]{fitdistr}} help 
document  carefully before setting the parameters in \code{fitdistr.args}.}

\item{wet.threshold}{The minimum value that is considered as a non-zero precipitation. Ignored when 
\code{precipitation = FALSE}. Default to 1 (assuming mm). See details on bias correction for precipitation.}

\item{n.quantiles}{Integer indicating the number of quantiles to be considered when method = "eqm", "dqm", "qdm". Default is NULL, 
that considers all quantiles, i.e. \code{n.quantiles = length(x[i,j])} (being \code{i} and \code{j} the 
coordinates in a single location).}

\item{extrapolation}{Character indicating the extrapolation method to be applied to correct values in
\code{newdata} that are out of the range of \code{x}. Extrapolation is applied only to the \code{"eqm"} method, 
thus, this argument is ignored if other bias correction method is selected. Default is \code{"none"} (do not extrapolate).}

\item{theta}{numeric indicating  upper threshold (and lower for the left tail of the distributions, if needed) 
above which precipitation (temperature) values are fitted to a Generalized Pareto Distribution (GPD). 
Values below this threshold are fitted to a gamma (normal) distribution. By default, 'theta' is the 95th 
percentile (and 5th percentile for the left tail). Only for \code{"gpqm"} method.}

\item{detrend}{logical. Detrend data prior to bias correction? Only for \code{"dqm"}. Default. TRUE.}

\item{isimip3.args}{Named list of arguments passed to function \code{\link{isimip3}}.}

\item{join.members}{Logical indicating whether members should be corrected independently (\code{FALSE}, the default),
or joined before performing the correction (\code{TRUE}). It applies to multimember grids only (otherwise ignored).}

\item{return.raw}{If TRUE, the nearest raw data to the observational reference is returned as the "var" dimension.
(Default to FALSE).}

\item{interpGrid.args}{Optional list fo the arguments passed to interpGrid. Configures the type of interpolation 
(Default "nearest") performed before bias adjustment.}

\item{parallel}{Logical. Should parallel execution be used?}

\item{max.ncores}{Integer. Upper bound for user-defined number of cores.}

\item{ncores}{Integer number of cores used in parallel computation. Self-selected number of
cores is used when \code{ncpus = NULL} (the default), or when \code{maxcores} exceeds the default \code{ncores} value.}
}
\value{
A calibrated grid of the same spatio-temporal extent than the input \code{"y"}
}
\description{
Implementation of several standard bias correction methods
}
\details{
The methods available are \code{"eqm"}, \code{"delta"}, 
\code{"scaling"}, \code{"pqm"}, \code{"gpqm"}, \code{"loci"}, 
\code{"ptr"}  (the four latter used only for precipitation), 
\code{"variance"} (only for temperature), \code{"dqm"} and \code{"qdm"}.

 These are next briefly described: 
 
\strong{Delta}

This method consists on adding to the observations the mean change signal (delta method).
This method is applicable to any kind of variable but it is preferable to avoid it for bounded variables
(e.g. precipitation, wind speed, etc.) because values out of the variable range could be obtained
(e.g. negative wind speeds...). This method corresponds to case g=1 and f=0 in Amengual et al. 2012. 

\strong{Scaling}

This method consists on scaling the simulation  with the difference (additive) or quotient (multiplicative) 
between the observed and simulated means in the train period. The \code{additive} or \code{multiplicative}
correction is defined by parameter \code{scaling.type} (default is \code{additive}).
The additive version is preferably applicable to unbounded variables (e.g. temperature) 
and the multiplicative to variables with a lower bound (e.g. precipitation, because it also preserves the frequency). 


\strong{eqm}

Empirical Quantile Mapping. This is a very extended bias correction method which consists on calibrating the simulated Cumulative Distribution Function (CDF) 
by adding to the observed quantiles both the mean delta change and the individual delta changes in the corresponding quantiles. 
This is equivalent to f=g=1 in Amengual et al. 2012. This method is applicable to any kind of variable.


\strong{pqm}

Parametric Quantile Mapping. It is based on the initial assumption that both observed and simulated intensity distributions are well approximated by a given distribution
(see \code{\link[MASS]{fitdistr}} to check available distributions), therefore is a parametric q-q map that uses the theorical instead of the empirical distribution.
For instance, the gamma distribution is described in Piani et al. 2010 and is applicable to precipitation. Other example is the weibull distribution, which
is applicable to correct wind data (Tie et al. 2014).
 
\strong{gpqm}
 
Generalized Quantile Mapping (described in Gutjahr and Heinemann 2013) is also a parametric quantile mapping (see
method 'pqm') but using two teorethical distributions, the gamma distribution and Generalized Pareto Distribution (GPD).
By default, It applies a gamma distribution to values under the threshold given by the 95th percentile 
(following Yang et al. 2010) and a general Pareto distribution (GPD) to values above the threshold. The threshold above 
which the GPD is fitted is the 95th percentile of the observed and the predicted wet-day distribution, respectively. If precip=FALSE
values below the 5th percentile of the observed and the predicted distributions are additionally fitted using GPD and 
the rest of the values of the distributions are fitted using a normal distribution.
The user can specify a different threshold(s) by modifying the parameter theta. 

\strong{mva}

Mean and Variance Adjustment.

\strong{variance}

Variance scaling of temperature. This method is described in Chen et al. 2011. It is applicable only to temperature. It corrects
the mean and variance of temperature time series.

\strong{loci}

Local intensity scaling of precipitation. This method is described in Schmidli et al. 2006. It adjust the mean as well as both wet-day frequencies and wet-day intensities.
The precipitation threshold is calculated such that the number of simulated days exceeding this threshold matches the number of observed days with precipitation larger than 1 mm.

\strong{ptr}

Power transformation of precipitation. This method is described in Leander and Buishand 2007 and is applicable only to precipitation. It adjusts the variance statistics of precipitation
time series in an exponential form. The power parameter is estimated on a monthly basis using a 90-day window centered on the interval. The power is defined by matching the coefficient
of variation of corrected daily simulated precipitation with the coefficient of variation of observed daily precipitation. It is calculated by root-finding algorithm using Brent's method.

\strong{dqm}

Detrended quantile matching with delta-method extrapolation, described in Cannon et al. 2015.
It consists on (i) removing the long-term mean (linear) trend; 
(ii) eqm is applied to the detrended series;
 (iii) the mean trend is then reapplied to the bias-adjusted series. 
It preserves the mean change signal in a climate change context.
It allows relative (multiplicative) and additive corrections.

\strong{qdm}

Quantile delta mapping, described in Cannon et al. 2015. 
It consists on (i) detrending the individual quantiles; 
(ii) QM is applied to the detrended series; 
(iii) the projected trends are then reapplied to the bias-adjusted quantiles.
It explicitly preserves the change signal in all quantiles. 
It allows relative (multiplicative) and additive corrections. 


\strong{isimip3}
}
\section{Note on the bias correction of precipitation}{


In the case of precipitation a frequency adaptation is performed in all versions of quantile mapping 
following Themeßl et al. (2012; https://doi.org/10.1007/s10584-011-0224-4), but sampling from the observed Gamma distribution instead of using 
linear interpolation. This is a preprocess to alleviate the problems arising when the dry day 
frequency in the raw model output is larger than in the observations. The opposite situation is 
automatically adjusted by quantile methods.
 
 The precipitation subroutines are switched-on when the variable name of the grid 
 (i.e., the value returned by \code{gridData$Variable$varName}) is one of the following: 
 \code{"pr"}, \code{"tp"} (this is the standard name defined in the vocabulary (\code{\link[cliamte4R.UDG]{C4R.vocabulary}}), \code{"precipitation"} or \code{"precip"}.
 Thus, caution must be taken to ensure that the correct bias correction is being undertaken when dealing with
 non-standard variables.
}

\examples{
\donttest{
require(climate4R.datasets)
data("EOBS_Iberia_pr")
data("CORDEX_Iberia_pr")
y <- EOBS_Iberia_pr
x <- CORDEX_Iberia_pr

# empirical
eqm1 <- biasCorrection(y = y, x = x,
                       precipitation = TRUE,
                       method = "eqm",
                       window = NULL,
                       wet.threshold = 0.1,
                       join.members = TRUE)
eqm1win <- biasCorrection(y = y, x = x,
                          precipitation = TRUE,
                          method = "eqm",
                          extrapolation = "none",
                          window = c(30, 15),
                          wet.threshold = 0.1)
eqm1folds <- biasCorrection(y = y, x = x,
                            precipitation = TRUE,
                            method = "eqm",
                            window = c(30, 15),
                            wet.threshold = 0.1,
                            cross.val = "kfold",
                            folds = list(1983:1989, 1990:1996, 1997:2002))


#parametric
pqm1.gamm <- biasCorrection(y = y, x = x,
                       method = "pqm",
                       precipitation = TRUE,
                       fitdistr.args = list(densfun = "gamma"))
pqm1.wei <- biasCorrection(y = y, x = x,
                       method = "pqm",
                       precipitation = TRUE,
                       fitdistr.args = list(densfun = "weibull"))

data("EOBS_Iberia_tas")
data("CORDEX_Iberia_tas")
y <- EOBS_Iberia_tas
x <- CORDEX_Iberia_tas
pqm1.norm <- biasCorrection(y = y, x = x,
           method = "pqm",
           fitdistr.args = list(densfun = "normal"))

# correction of future climate change data
data("CORDEX_Iberia_tas.rcp85")
newdata <- CORDEX_Iberia_tas.rcp85
eqm1win <- biasCorrection(y = y, x = x,
                          newdata = newdata,
                          method = "eqm",
                          extrapolation = "constant",
                          window = c(30, 15),
                          wet.threshold = 0.1)
pqm1.norm <- biasCorrection(y = y, x = x,
                       newdata = newdata,
                       method = "pqm",
                       fitdistr.args = list(densfun = "normal"))

# Correction of multimember datasets considering the joint
# distribution of all members
data("EOBS_Iberia_pr")
data("CFS_Iberia_pr")
y <- EOBS_Iberia_pr
x <- CFS_Iberia_pr
eqm.join <- biasCorrection(y = y, x = x,
                           precipitation = TRUE,
                           method = "eqm",
                           window = NULL,
                           wet.threshold = 0.1,
                           join.members = TRUE)
}
}
\references{
\itemize{
\item R.A.I. Wilcke, T. Mendlik and A. Gobiet (2013) Multi-variable error correction of regional climate models. Climatic Change, 120, 871-887

\item A. Amengual, V. Homar, R. Romero, S. Alonso, and C. Ramis (2012) A Statistical Adjustment of Regional Climate Model Outputs to Local Scales: Application to Platja de Palma, Spain. J. Clim., 25, 939-957

\item C. Piani, J. O. Haerter and E. Coppola (2009) Statistical bias correction for daily precipitation in regional climate models over Europe, Theoretical and Applied Climatology, 99, 187-192

\item O. Gutjahr and G. Heinemann (2013) Comparing precipitation bias correction methods for high-resolution regional climate simulations using COSMO-CLM, Theoretical and Applied Climatology, 114, 511-529

\item M. R. Tye, D. B. Stephenson, G. J. Holland and R. W. Katz (2014) A Weibull Approach for Improving Climate Model Projections of Tropical Cyclone Wind-Speed Distributions, Journal of Climate, 27, 6119-6133

\item Cannon, A.J., S.R. Sobie, and T.Q. Murdock (2015) Bias Correction of GCM Precipitation by Quantile Mapping: How Well Do Methods Preserve Changes in Quantiles and Extremes?. J. Climate, 28, 6938–6959, https://doi.org/10.1175/JCLI-D-14-00754.1
}
}
\seealso{
\code{\link{isimip}} for a trend-preserving method of model calibration.

Other downscaling: 
\code{\link{glimpr}()},
\code{\link{isimip}()}
}
\author{
S. Herrera, M. Iturbide, J. Bedia
}
\concept{downscaling}
