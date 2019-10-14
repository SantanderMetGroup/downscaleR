#     biasCorrection.R Bias correction methods
#
#     Copyright (C) 2017 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Bias correction methods
#' @description Implementation of several standard bias correction methods
#'
#' @template templateObsPredSim
#' @param method method applied. Current accepted values are \code{"eqm"}, \code{"delta"},
#'  \code{"scaling"}, \code{"pqm"} and \code{"gpqm"} \code{"variance"},\code{"loci"} and \code{"ptr"}. See details.
#' @param precipitation Logical for precipitation data (default to FALSE). If TRUE Adjusts precipitation 
#' frequency in 'x' (prediction) to the observed frequency in 'y'. This is a preprocess to bias correct 
#' precipitation data following Themeßl et al. (2012). To adjust the frequency, 
#' parameter \code{wet.threshold} is used (see below).
#' @param cross.val Logical (default to FALSE). Should cross-validation be performed? methods available are 
#' leave-one-out ("loo") and k-fold ("kfold") on an annual basis. The default option ("none") does not 
#' perform cross-validation.
#' @param folds Only requiered if \code{cross.val = "kfold"}. A list of vectors, each containing the years 
#' to be grouped in the corresponding fold.
#' @param wet.threshold The minimum value that is considered as a non-zero precipitation. Ignored when 
#' \code{precipitation = FALSE}. Default to 1 (assuming mm). See details on bias correction for precipitation.
#' @param window vector of length = 2 (or 1) specifying the time window width used to calibrate and the 
#' target days (days that are being corrected). The window is centered on the target day/s 
#' (window width >= target days). Default to \code{NULL}, which considers the whole period.
#' @param scaling.type Character indicating the type of the scaling method. Options are \code{"additive"} (default)
#' or \code{"multiplicative"} (see details). This argument is ignored if \code{"scaling"} is not 
#' selected as the bias correction method.
#' @param  fitdistr.args Further arguments passed to function \code{\link[MASS]{fitdistr}} 
#' (\code{densfun}, \code{start}, \code{...}). Only used when applying the "pqm" method 
#' (parametric quantile mapping). Please, read the \code{\link[MASS]{fitdistr}} help 
#' document  carefully before setting the parameters in \code{fitdistr.args}.
#' @param n.quantiles Integer indicating the number of quantiles to be considered when method = "eqm", "dqm", "qdm". Default is NULL, 
#' that considers all quantiles, i.e. \code{n.quantiles = length(x[i,j])} (being \code{i} and \code{j} the 
#' coordinates in a single location).
#' @param extrapolation Character indicating the extrapolation method to be applied to correct values in  
#' \code{newdata} that are out of the range of \code{x}. Extrapolation is applied only to the \code{"eqm"} method, 
#' thus, this argument is ignored if other bias correction method is selected. Default is \code{"none"} (do not extrapolate).
#' @param theta numeric indicating  upper threshold (and lower for the left tail of the distributions, if needed) 
#' above which precipitation (temperature) values are fitted to a Generalized Pareto Distribution (GPD). 
#' Values below this threshold are fitted to a gamma (normal) distribution. By default, 'theta' is the 95th 
#' percentile (5th percentile for the left tail). Only for \code{"gpqm"} method.
#' @param detrend logical. Detrend data prior to bias correction? Only for \code{"dqm"}. Default. TRUE.
#' @param join.members Logical indicating whether members should be corrected independently (\code{FALSE}, the default),
#'  or joined before performing the correction (\code{TRUE}). It applies to multimember grids only (otherwise ignored).
#' @template templateParallelParams
#'  
#' @details
#' 
#' The methods available are \code{"eqm"}, \code{"delta"}, 
#' \code{"scaling"}, \code{"pqm"}, \code{"gpqm"}, \code{"loci"}, 
#' \code{"ptr"}  (the four latter used only for precipitation), 
#' \code{"variance"} (only for temperature), \code{"dqm"} and \code{"qdm"}.
#' 
#'  These are next briefly described: 
#'  
#' \strong{Delta}
#'
#' This method consists on adding to the observations the mean change signal (delta method).
#' This method is applicable to any kind of variable but it is preferable to avoid it for bounded variables
#' (e.g. precipitation, wind speed, etc.) because values out of the variable range could be obtained
#' (e.g. negative wind speeds...). This method corresponds to case g=1 and f=0 in Amengual et al. 2012. 
#' 
#' \strong{Scaling}
#' 
#' This method consists on scaling the simulation  with the difference (additive) or quotient (multiplicative) 
#' between the observed and simulated means in the train period. The \code{additive} or \code{multiplicative}
#' correction is defined by parameter \code{scaling.type} (default is \code{additive}).
#' The additive version is preferably applicable to unbounded variables (e.g. temperature) 
#' and the multiplicative to variables with a lower bound (e.g. precipitation, because it also preserves the frequency). 
#' 
#' 
#' \strong{eqm}
#' 
#' Empirical Quantile Mapping. This is a very extended bias correction method which consists on calibrating the simulated Cumulative Distribution Function (CDF) 
#' by adding to the observed quantiles both the mean delta change and the individual delta changes in the corresponding quantiles. 
#' This is equivalent to f=g=1 in Amengual et al. 2012. This method is applicable to any kind of variable.
#' 
#' 
#' \strong{pqm}
#' 
#' Parametric Quantile Mapping. It is based on the initial assumption that both observed and simulated intensity distributions are well approximated by a given distribution
#' (see \code{\link[MASS]{fitdistr}} to check available distributions), therefore is a parametric q-q map that uses the theorical instead of the empirical distribution.
#' For instance, the gamma distribution is described in Piani et al. 2010 and is applicable to precipitation. Other example is the weibull distribution, which
#' is applicable to correct wind data (Tie et al. 2014).
#'  
#' \strong{gpqm}
#'  
#' Generalized Quantile Mapping (described in Gutjahr and Heinemann 2013) is also a parametric quantile mapping (see
#' method 'pqm') but using two teorethical distributions, the gamma distribution and Generalized Pareto Distribution (GPD).
#' By default, It applies a gamma distribution to values under the threshold given by the 95th percentile 
#' (following Yang et al. 2010) and a general Pareto distribution (GPD) to values above the threshold. the threshold above 
#' which the GPD is fitted is the 95th percentile of the observed and the predicted wet-day distribution, respectively. 
#' The user can specify a different threshold by modifying the parameter theta. It is applicable to precipitation data. 
#' 
#' \strong{mva}
#' 
#' Mean and Variance Adjustment.
#' 
#' \strong{variance}
#' 
#' Variance scaling of temperature. This method is described in Chen et al. 2011. It is applicable only to temperature. It corrects
#' the mean and variance of temperature time series.
#' 
#' \strong{loci}
#' 
#' Local intensity scaling of precipitation. This method is described in Schmidli et al. 2006. It adjust the mean as well as both wet-day frequencies and wet-day intensities.
#' The precipitation threshold is calculated such that the number of simulated days exceeding this threshold matches the number of observed days with precipitation larger than 1 mm.
#' 
#'\strong{ptr}
#'
#' Power transformation of precipitation. This method is described in Leander and Buishand 2007 and is applicable only to precipitation. It adjusts the variance statistics of precipitation
#' time series in an exponential form. The power parameter is estimated on a monthly basis using a 90-day window centered on the interval. The power is defined by matching the coefficient
#' of variation of corrected daily simulated precipitation with the coefficient of variation of observed daily precipitation. It is calculated by root-finding algorithm using Brent's method.
#'
#'\strong{dqm}
#'
#' Detrended quantile matching with delta-method extrapolation, described in Cannon et al. 2015.
#' It consists on (i) removing the long-term mean (linear) trend; 
#' (ii) eqm is applied to the detrended series;
#'  (iii) the mean trend is then reapplied to the bias-adjusted series. 
#' It preserves the mean change signal in a climate change context.
#' It allows relative (multiplicative) and additive corrections.
#'
#'\strong{qdm}
#'
#' Quantile delta mapping, described in Cannon et al. 2015. 
#' It consists on (i) detrending the individual quantiles; 
#' (ii) QM is applied to the detrended series; 
#' (iii) the projected trends are then reapplied to the bias-adjusted quantiles.
#' It explicitly preserves the change signal in all quantiles. 
#' It allows relative (multiplicative) and additive corrections. 
#'
#' @section Note on the bias correction of precipitation:
#' 
#' In the case of precipitation a frequency adaptation has been implemented in all versions of 
#' quantile mapping to alleviate the problems arising when the dry day frequency in the raw model output is larger
#'  than in the observations (Wilcke et al. 2013). 
#'  
#'  The precipitation subroutines are switched-on when the variable name of the grid 
#'  (i.e., the value returned by \code{gridData$Variable$varName}) is one of the following: 
#'  \code{"pr"}, \code{"tp"} (this is the standard name defined in the vocabulary (\code{\link[loadeR]{C4R.vocabulary}}), \code{"precipitation"} or \code{"precip"}.
#'  Thus, caution must be taken to ensure that the correct bias correction is being undertaken when dealing with
#'  non-standard variables.
#'     
#' 
#' @seealso \code{\link{isimip}} for a trend-preserving method of model calibration and \code{\link{quickDiagnostics}} 
#' for an outlook of the results.
#' @return A calibrated grid of the same spatio-temporal extent than the input \code{"y"}
#' @family downscaling
#' 
#' @importFrom transformeR redim subsetGrid getYearsAsINDEX getDim getWindowIndex
#' @importFrom abind adrop
#'
#' @references
#'
#' \itemize{
#' \item R.A.I. Wilcke, T. Mendlik and A. Gobiet (2013) Multi-variable error correction of regional climate models. Climatic Change, 120, 871-887
#'
#' \item A. Amengual, V. Homar, R. Romero, S. Alonso, and C. Ramis (2012) A Statistical Adjustment of Regional Climate Model Outputs to Local Scales: Application to Platja de Palma, Spain. J. Clim., 25, 939-957
#'
#' \item C. Piani, J. O. Haerter and E. Coppola (2009) Statistical bias correction for daily precipitation in regional climate models over Europe, Theoretical and Applied Climatology, 99, 187-192
#'
#' \item O. Gutjahr and G. Heinemann (2013) Comparing precipitation bias correction methods for high-resolution regional climate simulations using COSMO-CLM, Theoretical and Applied Climatology, 114, 511-529
#' 
#' \item M. R. Tye, D. B. Stephenson, G. J. Holland and R. W. Katz (2014) A Weibull Approach for Improving Climate Model Projections of Tropical Cyclone Wind-Speed Distributions, Journal of Climate, 27, 6119-6133
#' 
#' \item Cannon, A.J., S.R. Sobie, and T.Q. Murdock (2015) Bias Correction of GCM Precipitation by Quantile Mapping: How Well Do Methods Preserve Changes in Quantiles and Extremes?. J. Climate, 28, 6938–6959, https://doi.org/10.1175/JCLI-D-14-00754.1
#' }
#' @author S. Herrera, M. Iturbide, J. Bedia
#' @export
#' @examples \donttest{
#' data("EOBS_Iberia_pr")
#' data("CORDEX_Iberia_pr")
#' y <- EOBS_Iberia_pr
#' x <- CORDEX_Iberia_pr
#' 
#' # empirical
#' eqm1 <- biasCorrection(y = y, x = x,
#'                        precipitation = TRUE,
#'                        method = "eqm",
#'                        window = NULL,
#'                        wet.threshold = 0.1,
#'                        join.members = TRUE)
#' eqm1win <- biasCorrection(y = y, x = x,
#'                           precipitation = TRUE,
#'                           method = "eqm",
#'                           extrapolation = "none",
#'                           window = c(30, 15),
#'                           wet.threshold = 0.1)
#' eqm1folds <- biasCorrection(y = y, x = x,
#'                             precipitation = TRUE,
#'                             method = "eqm",
#'                             window = c(30, 15),
#'                             wet.threshold = 0.1,
#'                             cross.val = "kfold",
#'                             folds = list(1983:1989, 1990:1996, 1997:2002))
#' 
#' quickDiagnostics(y, x, eqm1, location = c(-2, 43))
#' quickDiagnostics(y, x, eqm1win, location = c(-2, 43))
#' quickDiagnostics(y, x, eqm1folds, location = c(-2, 43))
#' 
#' #parametric
#' pqm1.gamm <- biasCorrection(y = y, x = x,
#'                        method = "pqm",
#'                        precipitation = TRUE,
#'                        fitdistr.args = list(densfun = "gamma"))
#' quickDiagnostics(y, x, pqm1.gamm, location = c(-2, 43))
#' pqm1.wei <- biasCorrection(y = y, x = x,
#'                        method = "pqm",
#'                        precipitation = TRUE,
#'                        fitdistr.args = list(densfun = "weibull"))
#' quickDiagnostics(y, x, pqm1.wei, location = c(-2, 43))
#' 
#' data("EOBS_Iberia_tas")
#' data("CORDEX_Iberia_tas")
#' y <- EOBS_Iberia_tas
#' x <- CORDEX_Iberia_tas
#' pqm1.norm <- biasCorrection(y = y, x = x,
#'            method = "pqm",
#'            fitdistr.args = list(densfun = "normal"))
#' quickDiagnostics(y, x, pqm1.norm, location = c(-2, 43))
#' 
#' # correction of future climate change data
#' data("CORDEX_Iberia_tas.rcp85")
#' newdata <- CORDEX_Iberia_tas.rcp85
#' eqm1win <- biasCorrection(y = y, x = x,
#'                           newdata = newdata,
#'                           method = "eqm",
#'                           extrapolation = "constant",
#'                           window = c(30, 15),
#'                           wet.threshold = 0.1)
#' quickDiagnostics(y, x, eqm1win, location = c(-2, 43))
#' pqm1.norm <- biasCorrection(y = y, x = x,
#'                        newdata = newdata,
#'                        method = "pqm",
#'                        fitdistr.args = list(densfun = "normal"))
#' quickDiagnostics(y, x, pqm1.norm, location = c(-2, 43))
#' 
#' # Correction of multimember datasets considering the joint
#' # distribution of all members
#' data("EOBS_Iberia_pr")
#' data("CFS_Iberia_pr")
#' y <- EOBS_Iberia_pr
#' x <- CFS_Iberia_pr
#' eqm.join <- biasCorrection(y = y, x = x,
#'                            precipitation = TRUE,
#'                            method = "eqm",
#'                            window = NULL,
#'                            wet.threshold = 0.1,
#'                            join.members = TRUE)
#' }



biasCorrection <- function(y, x, newdata = NULL, precipitation = FALSE,
                           method = c("delta", "scaling", "eqm", "pqm", "gpqm", "loci","dqm","qdm"),
                           cross.val = c("none", "loo", "kfold"),
                           folds = NULL,
                           window = NULL,
                           scaling.type = c("additive", "multiplicative"),
                           fitdistr.args = list(densfun = "normal"),
                           wet.threshold = 1,
                           n.quantiles = NULL,
                           extrapolation = c("none", "constant"), 
                           theta = .95,
                           detrend=TRUE,
                           join.members = FALSE,
                           parallel = FALSE,
                           max.ncores = 16,
                           ncores = NULL) {
      if (method == "gqm") stop("'gqm' is not a valid choice anymore. Use method = 'pqm' instead and set fitdistr.args = list(densfun = 'gamma')")
      method <- match.arg(method, choices = c("delta", "scaling", "eqm", "pqm", "gpqm", "mva", "loci", "ptr", "variance","dqm","qdm"))
      cross.val <- match.arg(cross.val, choices = c("none", "loo", "kfold"))
      scaling.type <- match.arg(scaling.type, choices = c("additive", "multiplicative"))
      extrapolation <- match.arg(extrapolation, choices = c("none", "constant"))
      stopifnot(is.logical(join.members))
      nwdatamssg <- TRUE
      if (is.null(newdata)) {
            newdata <- x 
            nwdatamssg <- FALSE
      }
      if (cross.val == "none") {
            output <- biasCorrectionXD(y = y, x = x, newdata = newdata, 
                                       precipitation = precipitation,
                                       method = method,
                                       window = window,
                                       scaling.type = scaling.type,
                                       fitdistr.args = fitdistr.args,
                                       pr.threshold = wet.threshold, 
                                       n.quantiles = n.quantiles, 
                                       extrapolation = extrapolation, 
                                       theta = theta,
                                       join.members = join.members,
                                       detrend=detrend,
                                       parallel = parallel,
                                       max.ncores = max.ncores,
                                       ncores = ncores)
      } else {
            if (nwdatamssg) {
                  message("'newdata' will be ignored for cross-validation")
            }
            if (cross.val == "loo") {
                  years <- as.list(unique(getYearsAsINDEX(x)))
            } else if (cross.val == "kfold" & !is.null(folds)) {
                  years <- folds
            } else if (cross.val == "kfold" & is.null(folds)) {
                  stop("Fold specification is missing, with no default")
            }
            output.list <- lapply(1:length(years), function(i) {
                  target.year <- years[[i]]
                  rest.years <- setdiff(unlist(years), target.year)
                  station <- FALSE
                  if ("loc" %in% getDim(y)) station <- TRUE
                  yy <- redim(y, member = FALSE)
                  yy <- if (method == "delta") {
                        subsetGrid(yy, years = target.year, drop = FALSE)
                  } else {
                        subsetGrid(yy, years = rest.years, drop = FALSE)
                  }
                  if (isTRUE(station)) {
                        yy$Data <- adrop(yy$Data, drop = 3)
                        attr(yy$Data, "dimensions") <- c(setdiff(getDim(yy), c("lat", "lon")), "loc")
                  } else {
                        yy <- redim(yy, drop = TRUE)
                  }
                  newdata2 <- subsetGrid(x, years = target.year, drop = F)
                  xx <- subsetGrid(x, years = rest.years, drop = F)
                  message("Validation ", i, ", ", length(unique(years)) - i, " remaining")
                  biasCorrectionXD(y = yy, x = xx, newdata = newdata2, precipitation = precipitation,
                                   method = method,
                                   window = window,
                                   scaling.type = scaling.type,
                                   fitdistr.args = fitdistr.args,
                                   pr.threshold = wet.threshold, n.quantiles = n.quantiles, extrapolation = extrapolation, 
                                   theta = theta, join.members = join.members,
                                   detrend=detrend,
                                   parallel = parallel,
                                   max.ncores = max.ncores,
                                   ncores = ncores)
            })
            al <- which(getDim(x) == "time")
            Data <- sapply(output.list, function(n) unname(n$Data), simplify = FALSE)
            bindata <- unname(do.call("abind", c(Data, along = al)))
            output <- output.list[[1]]
            dimNames <- attr(output$Data, "dimensions")
            output$Data <- bindata
            attr(output$Data, "dimensions") <- dimNames
            output$Dates <- x$Dates
            output$Data[which(is.infinite(output$Data))] <- NA
      }
      return(output)
}


#' @keywords internal
#' @importFrom transformeR redim subsetGrid getDim getWindowIndex

biasCorrectionXD <- function(y, x, newdata, 
                             precipitation, 
                             method,
                             window,
                             scaling.type,
                             fitdistr.args,
                             pr.threshold, 
                             n.quantiles, 
                             extrapolation, 
                             theta,
                             join.members,
                             detrend,
                             parallel = FALSE,
                             max.ncores = 16,
                             ncores = NULL) {
      station <- FALSE
      if ("loc" %in% getDim(y)) station <- TRUE
      xy <- y$xyCoords
      suppressWarnings(suppressMessages(pred <- interpGrid(x, getGrid(y))))
      suppressWarnings(suppressMessages(sim <- interpGrid(newdata, getGrid(y))))
      delta.method <- method == "delta"
      precip <- precipitation
      message("[", Sys.time(), "] Argument precipitation is set as ", precip, ", please ensure that this matches your data.")
      bc <- y
      if (isTRUE(join.members) & getShape(redim(sim))["member"] > 1) {
            n.mem.aux <- getShape(sim)["member"]
            pred <- flatMemberDim(pred, station)
            pred <- redim(pred, drop = T)
            sim <- flatMemberDim(sim, station)
            sim <- redim(sim, drop = T)
            y <- bindGrid(rep(list(y), n.mem.aux), dimension = "time")
      } else if (isTRUE(join.members) & !getShape(redim(sim))["member"] > 1) {
            warning("There is only one member, argument join.members ignored.")
            join.members <- FALSE
      }
      y <- redim(y, drop = TRUE)
      y <- redim(y, member = FALSE, runtime = FALSE)
      pred <- redim(pred, member = TRUE, runtime = TRUE)
      sim <- redim(sim, member = TRUE, runtime = TRUE)
      dimNames <- attr(y$Data, "dimensions")
      n.run <- getShape(sim)["runtime"]
      n.mem <- getShape(sim)["member"]
      if (join.members & !is.null(window)) {
            message("[", Sys.time(), "] Window option is currently not supported for joined members and will be ignored")
            window <- NULL
      }
      if (!is.null(window)) {
            win <- getWindowIndex(y = y, x = pred, newdata = sim, window = window, delta.method = delta.method)
      } else {
            win <- list()
            indobservations <- match(as.POSIXct(pred$Dates$start), as.POSIXct(y$Dates$start))
            ## esto no mola, es para el caso especial de join members...hay que mirarlo
            if (length(indobservations) > length(unique(indobservations))) indobservations <- 1:length(indobservations) 
            win[["Window1"]] <- list("obsWindow" = indobservations, "window" = 1:getShape(pred)["time"], "step" = 1:getShape(sim)["time"])
            if (delta.method) win[["Window1"]][["deltaind"]] <- indobservations
      }
      message("[", Sys.time(), "] Number of windows considered: ", length(win), "...")
      winarr <- array(dim = dim(sim$Data))
      if (delta.method) winarr <- array(dim = c(n.run, n.mem, getShape(y)))
      for (j in 1:length(win)) {
            yind <- win[[j]]$obsWindow
            outind <- win[[j]]$step
            if (delta.method) {
                  yind <- win[[j]]$deltaind
                  outind <- win[[j]]$deltaind
            } 
            yw <- y$Data[yind,,, drop = FALSE]
            pw <- pred$Data[,,win[[j]]$window,,, drop = FALSE]
            sw <- sim$Data[,,win[[j]]$step,,, drop = FALSE]
            runarr <- lapply(1:n.run, function(l){
                  memarr <- lapply(1:n.mem, function(m){
                        #join members message
                        if (j == 1 & m == 1) {
                              if (!isTRUE(join.members)) {
                                    message("[", Sys.time(), "] Bias-correcting ", n.mem, " members separately...")
                              } else {
                                    message("[", Sys.time(), "] Bias-correcting ", attr(pred, "orig.mem.shape"), " members considering their joint distribution...")
                              }
                        }
                        o = yw[, , , drop = FALSE]
                        p = adrop(pw[l, m, , , , drop = FALSE], drop = c(T, T, F, F, F))
                        s = adrop(sw[l, m, , , , drop = FALSE], drop = c(T, T, F, F, F))
                        data <- list(o, p, s)
                        if (!station) {
                              data <- lapply(1:length(data), function(x) {
                                    attr(data[[x]], "dimensions") <- dimNames
                                    abind(array3Dto2Dmat(data[[x]]), along = 3)
                              }) 
                        }
                        o <- lapply(seq_len(ncol(data[[1]])), function(i) data[[1]][,i,1])
                        p <- lapply(seq_len(ncol(data[[2]])), function(i) data[[2]][,i,1])
                        s <- lapply(seq_len(ncol(data[[3]])), function(i) data[[3]][,i,1])
                        mat <- biasCorrection1D(o, p, s,
                                                method = method,
                                                scaling.type = scaling.type,
                                                fitdistr.args = fitdistr.args,
                                                precip = precip,
                                                pr.threshold = pr.threshold,
                                                n.quantiles = n.quantiles,
                                                extrapolation = extrapolation,
                                                theta = theta,
                                                detrend=detrend,
                                                parallel = parallel,
                                                max.ncores = max.ncores,
                                                ncores = ncores)  
                        if (!station) mat <- mat2Dto3Darray(mat, xy$x, xy$y)
                        mat
                  })
                  unname(do.call("abind", list(memarr, along = 0)))
            })
            winarr[,,outind,,] <- unname(do.call("abind", list(runarr, along = 0))) 
      }
      bc$Data <- unname(do.call("abind", list(winarr, along = 3)))
      attr(bc$Data, "dimensions") <- attr(sim$Data, "dimensions")
      if (station) bc <- redim(bc, loc = TRUE)
      bc$Dates <- sim$Dates
      ## Recover the member dimension when join.members=TRUE:
      if (isTRUE(join.members)) {
            if (method == "delta") {
                  bc <- recoverMemberDim(plain.grid = pred, bc.grid = bc, newdata = newdata)
            }else{
                  bc <- recoverMemberDim(plain.grid = sim, bc.grid = bc, newdata = newdata)      
            }
      } else {
            bc$InitializationDates <- sim$InitializationDates
            bc$Members <- sim$Members
      }
      attr(bc$Variable, "correction") <- method
      bc <- redim(bc, drop = TRUE)
      message("[", Sys.time(), "] Done.")
      return(bc)
}

#' @title Bias correction methods on 1D data
#' @description Implementation of several standard bias correction methods
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @param method method applied. Current accepted values are \code{"eqm"}, \code{"delta"},
#'  \code{"scaling"}, \code{"pqm"} , \code{"gpqm"}, \code{"mva"}, \code{"variance"},\code{"loci"} and \code{"ptr"}. See details in 
#'  function \code{\link{biasCorrection}}.
#' @param scaling.type Character indicating the type of the scaling method. Options are \code{"additive"} 
#' or \code{"multiplicative"} (see details). This argument is ignored if \code{"scaling"} is not 
#' selected as the bias correction method.
#' @param  fitdistr.args Further arguments passed to function \code{\link[MASS]{fitdistr}} 
#' (\code{densfun}, \code{start}, \code{...}). Only used when applying the "pqm" method 
#' (parametric quantile mapping). Please, read the \code{\link[MASS]{fitdistr}} help 
#' document  carefully before setting the parameters in \code{fitdistr.args}.
#' @param precip Logical for precipitation data. If TRUE Adjusts precipitation 
#' frequency in 'x' (prediction) to the observed frequency in 'y'. This is a preprocess to bias correct 
#' precipitation data following Themeßl et al. (2012). To adjust the frequency, 
#' parameter \code{pr.threshold} is used (see below).
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored when 
#' \code{precip = FALSE}. See details in function \code{biasCorrection}.
#' @param n.quantiles Integer indicating the number of quantiles to be considered when method = "eqm". 
#' @param extrapolation Character indicating the extrapolation method to be applied to correct values in  
#' \code{newdata} that are out of the range of \code{x}. Extrapolation is applied only to the \code{"eqm"} method, 
#' thus, this argument is ignored if other bias correction method is selected.
#' @param theta numeric indicating  upper threshold (and lower for the left tail of the distributions, if needed) 
#' above which precipitation (temperature) values are fitted to a Generalized Pareto Distribution (GPD). 
#' Values below this threshold are fitted to a gamma (normal) distribution. By default, 'theta' is the 95th 
#' percentile (5th percentile for the left tail). Only for \code{"gpqm"} method.
#' @param detrend logical. Detrend data prior to bias correction? Only for \code{"dqm"}. Default. TRUE.
#' @template templateParallelParams
#' 
#' 
#' @importFrom transformeR parallelCheck selectPar.pplyFun
#' @keywords internal
#' @author M. Iturbide

biasCorrection1D <- function(o, p, s,
                             method, 
                             scaling.type,
                             fitdistr.args,
                             precip, 
                             pr.threshold,
                             n.quantiles,
                             extrapolation, 
                             theta,
                             detrend,
                             parallel = FALSE,
                             max.ncores = 16,
                             ncores = NULL) {
      parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
      mapply_fun <- selectPar.pplyFun(parallel.pars, .pplyFUN = "mapply")
      if (parallel.pars$hasparallel) on.exit(parallel::stopCluster(parallel.pars$cl))
      if (method == "delta") {
            mapply_fun(delta, o, p, s)
      } else if (method == "scaling") {
            mapply_fun(scaling, o, p, s, MoreArgs = list(scaling.type = scaling.type))
      } else if (method == "eqm") {
            suppressWarnings(
                  mapply_fun(eqm, o, p, s, MoreArgs = list(precip, pr.threshold, n.quantiles, extrapolation))
            )
      } else if (method == "pqm") {
            suppressWarnings(
            mapply_fun(pqm, o, p, s, MoreArgs = list(fitdistr.args, precip, pr.threshold))
            )
      } else if (method == "gpqm") {
            mapply_fun(gpqm, o, p, s, MoreArgs = list(precip, pr.threshold, theta))
      } else if (method == "mva") {
            mapply_fun(mva, o, p, s) 
      } else if (method == "variance") {
            mapply_fun(variance, o, p, s, MoreArgs = list(precip))
      } else if (method == "loci") {
            mapply_fun(loci, o, p, s, MoreArgs = list(precip, pr.threshold))
      } else if (method == "ptr") {
            mapply_fun(ptr, o, p, s, MoreArgs = list(precip))
      } else if (method == "dqm") {
            mapply_fun(dqm, o, p, s, MoreArgs = list(precip, pr.threshold, n.quantiles, detrend))
      } else if (method == "qdm") {
            mapply_fun(qdm, o, p, s, MoreArgs = list(precip, pr.threshold, n.quantiles))
      }
      #INCLUIR AQUI METODOS NUEVOS
}

#' @title adjustPrecipFreq
#' @description Adjusts precipitation frequency in 'p' (prediction) to the observed frequency in 'o'. 
#' It constitutes a preprocess to bias correct precipitation data following Themeßl et al. (2012). 
#' @param obs A vector (e.g. station data) containing the observed climate data for the training period
#' @param pred A vector containing the simulated climate by the model for the training period. 
#' @param threshold The minimum value that is considered as a non-zero precipitation. 
#' @importFrom MASS fitdistr
#' @keywords internal
#' @importFrom stats rgamma
#' @author S. Herrera and M. Iturbide

adjustPrecipFreq <- function(obs, pred, threshold){
      o <- obs[!is.na(obs)]
      p <- pred[!is.na(pred)]
      # Number of dry days in 'o' 
      nPo <- sum(as.double(o < threshold))
      # Number of dry days that must be in 'p' to equal precip frequency in 'o'
      nPp <- ceiling(length(p) * nPo / length(o))
      # Index and values of ordered 'p'
      ix <- sort(p, decreasing = FALSE, index.return = TRUE)$ix
      Ps <- sort(p, decreasing = FALSE)
      Pth <- max(Ps[nPp:(nPp + 1)], na.rm = TRUE) # in case nPp == length(Ps)
      # Themeßl (Themessl) modification (simulating rain for model dry days) 
      inddrzl <- which(Ps[(nPp + 1):length(Ps)] < threshold)
      if (length(inddrzl) > 0) { 
            Os <- sort(o, decreasing = FALSE, na.last = NA)
            indO <- ceiling(length(Os) * (nPp + max(inddrzl))/length(Ps))
            auxOs <- Os[(nPo + 1):indO]
            if (length(unique(auxOs)) > 6) {
                  # simulate precip for 'p' with a gamma adjusted in 'o' for values between
                  auxGamma <- fitdistr(auxOs, "gamma")
                  Ps[(nPp + 1):(nPp + max(inddrzl))] <- rgamma(length(inddrzl), auxGamma$estimate[1], rate = auxGamma$estimate[2])
            } else {
                  Ps[(nPp + 1):(nPp + max(inddrzl))] <- mean(auxOs)
            }
            # order 'Ps' after simulation
            Ps <- sort(Ps, decreasing = FALSE, na.last = NA)
      }
      # Make 0-s
      if (nPo > 0) {
            ind <- min(nPp, length(p))
            Ps[1:ind] <- 0
      }
      p[ix] <- Ps
      pred[!is.na(pred)] <- p
      return(list("nP" = c(nPo,nPp), "Pth" = Pth, "p" = pred)) 
}
#end

#' @title Delta method for bias correction
#' @description Implementation of Delta method for bias correction 
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{x}, but considering the test period.
#' @keywords internal
#' @author S. Herrera and M. Iturbide

delta <- function(o, p, s){
      corrected <- o + (mean(s) - mean(p))
      return(corrected)
}


#' @title Scaling method for bias correction
#' @description Implementation of Scaling method for bias correction 
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{x}, but considering the test period.
#' @param scaling.type Character indicating the type of the scaling method. Options are \code{"additive"} (default)
#' or \code{"multiplicative"} (see details). This argument is ignored if \code{"scaling"} is not selected as the bias correction method.
#' @keywords internal
#' @author S. Herrera and M. Iturbide

scaling <- function(o, p, s, scaling.type){
      if (scaling.type == "additive") {
            s - mean(p) + mean(o, na.rm = TRUE)
      } else if (scaling.type == "multiplicative") {
            (s/mean(p)) * mean(o, na.rm = TRUE)
      }
}


#' @title Parametric Quantile Mapping method for bias correction
#' @description Implementation of Parametric Quantile Mapping method for bias correction 
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{x}, but considering the test period.
#' @param  fitdistr.args Further arguments passed to function \code{\link[MASS]{fitdistr}} 
#' (\code{densfun}, \code{start}, \code{...}). Only used when applying the "pqm" method 
#' (parametric quantile mapping). Please, read the \code{\link[MASS]{fitdistr}} help 
#' document  carefully before parameter setting in \code{fitdistr.args}.
#' @param precip Logical for precipitation data. If TRUE Adjusts precipitation 
#' frequency in 'x' (prediction) to the observed frequency in 'y'. This is a preprocess to bias correct 
#' precipitation data following Themeßl et al. (2012). To adjust the frequency, 
#' parameter \code{pr.threshold} is used (see below).
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored when 
#' \code{precip = FALSE}. See details in function \code{biasCorrection}.
#' @importFrom MASS fitdistr
#' @importFrom stats pgamma qgamma 
#' @keywords internal
#' @author S. Herrera and M. Iturbide

pqm <- function(o, p, s, fitdistr.args, precip, pr.threshold){
      dfdistr <- cbind("df" = c( "beta", "cauchy", "chi-squared", "exponential", "f", "gamma", "geometric", "log-normal", "lognormal", "logistic", "negative binomial", "normal", "Poisson", "t", "weibull"),
                        "p" = c("pbeta", "pcauchy", "pchisq", "pexp", "pf", "pgamma", "pegeom", "plnorm", "plnorm", "plogis", "pnbinom", "pnorm", "ppois", "pt", "pweibull"),
                        "q" = c("qbeta", "qcauchy", "qchisq", "qexp", "qf", "qgamma", "qegeom", "qlnorm", "qlnorm", "qlogis", "qnbinom", "qnorm", "qpois", "qt", "qweibull"))
      fitdistr.args <- fitdistr.args[which(names(fitdistr.args) != "x")]
      statsfunp <- unname(dfdistr[which(dfdistr[,"df"] == fitdistr.args$densfun), "p"])
      statsfunq <- unname(dfdistr[which(dfdistr[,"df"] == fitdistr.args$densfun), "q"])
      run <- TRUE
      ind.o <- 1:length(o)
      ind.p <- 1:length(p)
      rain <- 1:length(s)
      if (precip) {
            threshold <- pr.threshold
            if (any(!is.na(o))) {
                  params <-  adjustPrecipFreq(o, p, threshold)
                  p <- params$p
                  nP <- params$nP
                  Pth <- params$Pth
            } else {
                  nP = NULL
            }
            if (is.null(nP)) {
                  run <- FALSE
                  s <- rep(NA, length(s))
            } else if (nP[1] < length(o)) {
                  ind.o <- which(o > threshold & !is.na(o))
                  ind.p <- which(p > 0 & !is.na(p))
                  rain <- which(s > Pth & !is.na(s))
                  noRain <- which(s <= Pth & !is.na(s))
            } else {
                  run <- FALSE
                  warning("For the window step selected, location without rainfall above the threshold.\n no bias correction applied in location.")
            } 
      }
      if (all(is.na(o[ind.o]))) {
            run <- FALSE
            s <- rep(NA, length(s))
      }
      if (run) {
            fitdistr.args.o <- c("x" = list(o[ind.o]), fitdistr.args)
            fitdistr.args.p <- c("x" = list(p[ind.p]), fitdistr.args)
            obsGamma <- tryCatch({do.call("fitdistr", fitdistr.args.o)}, error = function(err){NULL})
            prdGamma <- tryCatch({do.call("fitdistr", fitdistr.args.p)}, error = function(err){NULL})
            if (!is.null(prdGamma) & !is.null(obsGamma)) {
                  statsfun.args <- c(list(s[rain]), as.list(prdGamma$estimate))
                  auxF <- do.call(statsfunp, statsfun.args)
                  statsfun.args <- c(list(auxF), as.list(obsGamma$estimate))
                  s[rain] <- do.call(statsfunq, statsfun.args)
                  if (precip) s[noRain] <- 0
            } else {
                  warning("Fitting error for location and selected 'densfun'.")
            }
      }   
      return(s)      
}   
#end

#' @title Empirical Quantile Mapping method for bias correction
#' @description Implementation of Empirical Quantile Mapping method for bias correction 
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @param precip Logical for precipitation data. If TRUE Adjusts precipitation 
#' frequency in 'x' (prediction) to the observed frequency in 'y'. This is a preprocess to bias correct 
#' precipitation data following Themeßl et al. (2012). To adjust the frequency, 
#' parameter \code{pr.threshold} is used (see below).
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored when 
#' \code{precip = FALSE}. See details in function \code{biasCorrection}.
#' @param n.quantiles Integer indicating the number of quantiles to be considered when method = "eqm". Default is NULL, 
#' that considers all quantiles, i.e. \code{n.quantiles = length(p)}.
#' @param extrapolation Character indicating the extrapolation method to be applied to correct values in  
#' \code{"s"} that are out of the range of \code{"p"}. Extrapolation is applied only to the \code{"eqm"} method, 
#' thus, this argument is ignored if other bias correction method is selected. 
#' @keywords internal
#' @importFrom stats approxfun ecdf quantile
#' @author S. Herrera and M. Iturbide

eqm <- function(o, p, s, precip, pr.threshold, n.quantiles, extrapolation){
      if (precip == TRUE) {
            threshold <- pr.threshold
            if (any(!is.na(o))) {
                  params <-  adjustPrecipFreq(o, p, threshold)
                  p <- params$p
                  nP <- params$nP
                  Pth <- params$Pth
            } else {
                  nP = NULL
            }
            smap <- rep(NA, length(s))
            if (any(!is.na(p)) & any(!is.na(o))) {
                  if (length(which(p > Pth)) > 0) { 
                        noRain <- which(s <= Pth & !is.na(s))
                        rain <- which(s > Pth & !is.na(s))
                        drizzle <- which(s > Pth  & s  <= min(p[which(p > Pth)], na.rm = TRUE) & !is.na(s))
                        if (length(rain) > 0) {
                              eFrc <- tryCatch({ecdf(s[rain])}, error = function(err) {stop("There are not precipitation days in newdata for the step length selected in one or more locations. Try to enlarge the window step")})
                              if (is.null(n.quantiles)) n.quantiles <- length(p)
                              bins <- n.quantiles
                              qo <- quantile(o[which(o > threshold & !is.na(o))], prob = seq(1/bins,1 - 1/bins,1/bins), na.rm = T)
                              qp <- quantile(p[which(p > Pth)], prob = seq(1/bins,1 - 1/bins,1/bins), na.rm = T)
                              p2o <- tryCatch({approxfun(qp, qo, method = "linear")}, error = function(err) {NA})
                              smap <- s
                              smap[rain] <- if (suppressWarnings(!is.na(p2o))) {
                                    p2o(s[rain])
                              }else{
                                    s[rain] <- NA
                              }
                              # Linear extrapolation was discarded due to lack of robustness 
                              if (extrapolation == "constant") {
                                    smap[rain][which(s[rain] > max(qp, na.rm = TRUE))] <- s[rain][which(s[rain] > max(qp, na.rm = TRUE))] + (qo[length(qo)] - qp[length(qo)])
                                    smap[rain][which(s[rain] < min(qp, na.rm = TRUE))] <- s[rain][which(s[rain] < min(qp, na.rm = TRUE))] + (qo[1] - qp[1]) 
                              } else {
                                    smap[rain][which(s[rain] > max(qp, na.rm = TRUE))] <- qo[length(qo)]
                                    smap[rain][which(s[rain] < min(qp, na.rm = TRUE))] <- qo[1]
                              }
                        }else{
                              smap <- rep(0, length(s))
                              warning("There are not precipitation days in newdata for the step length selected in one or more locations. Consider the possibility of enlarging the window step")
                        }
                        if (length(drizzle) > 0) {
                              smap[drizzle] <- quantile(s[which(s > min(p[which(p > Pth)], na.rm = TRUE) & !is.na(s))], probs = eFrc(s[drizzle]), na.rm = TRUE, type = 4)
                        }
                        smap[noRain] <- 0
                  } else { ## For dry series
                        smap <- s
                        warning('No rainy days in the prediction. Bias correction is not applied') 
                  }
            }
      } else {
            if (all(is.na(o))) {
                  smap <- rep(NA, length(s))
            } else if (all(is.na(p))) {
                  smap <- rep(NA, length(s))
            }else if (any(!is.na(p)) & any(!is.na(o))) {
                  if (is.null(n.quantiles)) n.quantiles <- length(p)
                  bins <- n.quantiles
                  qo <- quantile(o, prob = seq(1/bins,1 - 1/bins,1/bins), na.rm = TRUE)
                  qp <- quantile(p, prob = seq(1/bins,1 - 1/bins,1/bins), na.rm = TRUE)
                  p2o <- approxfun(qp, qo, method = "linear")
                  smap <- p2o(s)
                  if (extrapolation == "constant") {
                        smap[which(s > max(qp, na.rm = TRUE))] <- s[which(s > max(qp, na.rm = TRUE))] + (qo[length(qo)] - qp[length(qo)])
                        smap[which(s < min(qp, na.rm = TRUE))] <- s[which(s < min(qp, na.rm = TRUE))] + (qo[1] - qp[1]) 
                  } else {
                        smap[which(s > max(qp, na.rm = TRUE))] <- qo[length(qo)]
                        smap[which(s < min(qp, na.rm = TRUE))] <- qo[1]
                  }
            } 
      }
      return(smap)
}
#end

#' @title Generalized Quantile Mapping method for bias correction
#' @description Implementation of Generalized Quantile Mapping method for bias correction 
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @param precip Logical for precipitation data. If TRUE Adjusts precipitation 
#' frequency in 'x' (prediction) to the observed frequency in 'y'. This is a preprocess to bias correct 
#' precipitation data following Themeßl et al. (2012). To adjust the frequency, 
#' parameter \code{pr.threshold} is used (see below).
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored when 
#' \code{precip = FALSE}. See details in function \code{biasCorrection}.
#' @param theta numeric indicating  upper threshold (and lower for the left tail of the distributions, if needed) 
#' above which precipitation (temperature) values are fitted to a Generalized Pareto Distribution (GPD). 
#' Values below this threshold are fitted to a gamma (normal) distribution. By default, 'theta' is the 95th 
#' percentile (5th percentile for the left tail). 
#' @importFrom evd fpot
#' @importFrom MASS fitdistr
#' @importFrom evd qgpd pgpd
#' @importFrom stats quantile pgamma qgamma
#' @keywords internal
#' @author S. Herrera and M. Iturbide

gpqm <- function(o, p, s, precip, pr.threshold, theta) { 
      if (precip == FALSE) {
            stop("method gpqm is only applied to precipitation data")
      } else {
            threshold <- pr.threshold
            if (any(!is.na(o))) {
                  params <-  adjustPrecipFreq(o, p, threshold)
                  p <- params$p
                  nP <- params$nP
                  Pth <- params$Pth
            } else {
                  nP = NULL
            }
            if (is.null(nP)) {
                  s <- rep(NA, length(s))
            } else if (nP[1] < length(o)) {
                  ind <- which(o > threshold & !is.na(o))
                  indgamma <- ind[which(o[ind] < quantile(o[ind], theta))]
                  indpareto <- ind[which(o[ind] >= quantile(o[ind], theta))]
                  obsGQM <- fitdistr(o[indgamma],"gamma")
                  obsGQM2 <- fpot(o[indpareto], quantile(o[ind], theta), "gpd", std.err = FALSE)
                  ind <- which(p > 0 & !is.na(p))
                  indgammap <- ind[which(p[ind] < quantile(p[ind],theta))]
                  indparetop <- ind[which(p[ind] >= quantile(p[ind], theta))]
                  prdGQM <- fitdistr(p[indgammap], "gamma")
                  prdGQM2 <- fpot(p[indparetop], quantile(p[ind], theta), "gpd", std.err = FALSE)
                  rain <- which(s > Pth & !is.na(s))
                  noRain <- which(s <= Pth & !is.na(s))
                  indgammasim <- rain[which(s[rain] < quantile(p[ind], theta))]
                  indparetosim <- rain[which(s[rain] >= quantile(p[ind], theta))]
                  auxF <- pgamma(s[indgammasim], prdGQM$estimate[1], rate = prdGQM$estimate[2])
                  auxF2 <- pgpd(s[indparetosim], loc = 0, scale = prdGQM2$estimate[1], shape = prdGQM2$estimate[2])
                  s[indgammasim] <- qgamma(auxF, obsGQM$estimate[1], rate = obsGQM$estimate[2])
                  s[indparetosim[which(auxF2 < 1)]] <- qgpd(auxF2[which(auxF2 < 1)], loc = 0, scale = obsGQM2$estimate[1], shape = obsGQM2$estimate[2])
                  s[indparetosim[which(auxF2 == 1)]] <- max(o[indpareto], na.rm = TRUE)
                  s[noRain] <- 0
            } else {
                  warning("There is at least one location without rainfall above the threshold.\n In this (these) location(s) none bias correction has been applied.")
            }  
      }
      return(s)
}

#end

#' @title Mean and Variance Adjustment
#' @description Mean and Variance Adjustment method for bias correction
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @keywords internal
#' @author M. Iturbide

mva <- function(o, p, s){
      corrected <- (s - mean(p)) + sd(o)/sd(p) + mean(o)
      return(corrected)
}


#' @title Variance scaling of temperature
#' @description Implementation of Variance scaling of temperature method for bias correction
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @param precip Logical indicating if o, p, s is temperature data.
#' @keywords internal
#' @author B. Szabo-Takacs

variance <- function(o, p, s, precip) {
      if (precip == FALSE) {
            t_dif <- mean(o, na.rm = TRUE) - mean(p, na.rm = TRUE)
            t1 <- p + rep(t_dif, length(p), 1)
            t1_m <- mean(t1,na.rm = TRUE) 
            t2 <- t1 - rep(t1_m,length(t1),1)
            o_s <- sd(o,na.rm = TRUE) 
            t2_s <- sd(t2,na.rm = TRUE) 
            tsig <- o_s/t2_s
            t1 <- t1_m <- t2 <- o_s <- t2_s <- NULL
            t1 <- s + rep(t_dif, length(s), 1)
            t1_m <- mean(t1, na.rm = TRUE)
            t2 <- t1 - rep(t1_m, length(t1), 1)
            t3 <- t2 * rep(tsig, length(t2), 1)
            tC <- t3 + rep(t1_m, length(t3), 1)
            t1 <- t1_m <- t2 <- t3 <- NULL
            return(tC)
      } else {
            stop("method variance is only applied to temperature data")
      }
}


#' @title Local intensity scaling of precipitation
#' @description Implementation of Local intensity scaling of precipitation method for bias correction based on Vincent Moron's local_scaling function in weaclim toolbox in Matlab
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training or test period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @param precip Logical indicating if o, p, s is precipitation data.
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored when 
#' \code{precip = FALSE}. See details in function \code{biasCorrection}.
#' @author B. Szabo-Takacs

loci <- function(o, p, s, precip, pr.threshold){
      if (precip == FALSE) { 
            stop("method loci is only applied to precipitation data")
      } else {
            threshold <- pr.threshold
            l <- length(which(o > threshold))
            gcmr <- rev(sort(p))
            gcmrs <- rev(sort(s))
            Pgcm <- gcmr[l + 1]
            Pgcms <- gcmrs[l + 1]
            mobs <- mean(o[which(o > threshold)], na.rm = TRUE)
            mgcm <- mean(p[which(p > Pgcm)], na.rm = TRUE)
            scaling <- (mobs - threshold) / (mgcm - Pgcm)
            GCM <- (scaling*(s - Pgcms)) + threshold
            GCM[which(GCM < threshold)] <- 0
      }
      return(GCM)
}


#' @title Power transformation of precipitation
#' @description Implementation of Power transformation of precipitation method for bias correction 
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @param precip Logical indicating if o, p, s is precipitation data.
#' @importFrom stats uniroot
#' @keywords internal
#' @author S. Herrera and B. Szabo-Takacs

ptr <- function(o, p, s, precip) {
      if (precip == FALSE) { 
            stop("method power transformation is only applied to precipitation data")
      } else {
            b <- NaN
            cvO <- sd(o,na.rm = TRUE) / mean(o, na.rm = TRUE)
            if (!is.na(cvO)) {
                  bi <- try(uniroot(function(x)
                        varCoeficient(x, abs(p), cvO), c(0,1), extendInt = "yes"), silent = TRUE)
                  if ("try-error" %in% class(bi)) {  # an error occurred
                        b <- NA
                  } else {
                        b <- bi$root
                  }
            }
            p[p < 0] <-  0
            s[s < 0] <-  0
            aux_c <- p^rep(b,length(p),1)
            aux <- s^rep(b,length(s),1)
            prC <- aux * rep((mean(o, na.rm = TRUE) / mean(aux_c, na.rm = TRUE)), length(s), 1)
            aux <- aux_c <- NULL
      }
      return(prC)
}


#' @title VarCoeficient
#' @description preprocess to power transformation of precipitation
#' @param delta A vector of power parameter
#' @param data A vector containing the simulated climate by the model for training period
#' @param cv A vector containing coefficient of variation of observed climate data
#' @keywords internal
#' @author S. Herrera and B. Szabo-Takacs

varCoeficient <- function(delta,data,cv){
      y <- cv - sd((data^delta), na.rm = TRUE)/mean((data^delta), na.rm = TRUE)
      return(y)
}


#' @title Concatenate members
#' @description Concatenate members as a single time series for using their joint distribution in bias correction
#' @param grid Input (multimember) grid
#' @return A grid without members, with additional attributes to retrieve the original structure after bias correction
#' @seealso \code{\link{recoverMemberDim}}, for recovering the original structure after bias correction.
#' @keywords internal
#' @importFrom transformeR subsetGrid redim getShape bindGrid
#' @author J Bedia

flatMemberDim <- function(grid, station) {
      grid <- redim(grid, member = TRUE, loc = station)     
      n.mem.join <- getShape(grid, "member")
      n.time.join <- getShape(grid, "time")
      aux.ltime <- lapply(1:n.mem.join, function(x) {
            subsetGrid(grid, members = x)
      })
      out <- do.call("bindGrid", c(aux.ltime, dimension = "time"))
      attr(out, "orig.mem.shape") <- n.mem.join
      attr(out, "orig.time.shape") <- n.time.join
      return(out)
}

#' @title Recover member multimember structure
#' @description Recover member multimember structure after application of \code{\link{flatMemberDim}}
#' @param plain.grid A \dQuote{flattened} grid used as predictor in \code{biasCorrection} (the 'pred' object)
#' @param bc.grid The bias-corrected output (the 'bc' object), still without its member structure 
#' @param newdata The 'newdata' object, needed to recover relevant metadata (i.e. initialization dates and member names)
#' @return A (bias-corrected) multimember grid
#' @keywords internal
#' @importFrom transformeR subsetDimension bindGrid
#' @seealso \code{\link{flatMemberDim}}, for \dQuote{flattening} the member structure
#' @author J Bedia

recoverMemberDim <- function(plain.grid, bc.grid, newdata) {
      bc <- bc.grid
      nmem <- attr(plain.grid, "orig.mem.shape")
      ntimes <- attr(plain.grid, "orig.time.shape")
      # bc$Dates <- lapply(bc$Dates, "rep", nmem)
      aux.list <- lapply(1:nmem, function(m) {
            aux <- subsetDimension(grid = bc, dimension = "time", indices = seq(m, nmem * ntimes, nmem))
            aux$InitializationDates <- newdata$InitializationDates[[m]]
            aux$Members <- newdata$Members[[m]]
            return(aux)
      })
      do.call("bindGrid", c(aux.list, dimension = "member"))
}

#' @title Detrended quantile matching 
#' @description Detrended quantile matching with delta-method extrapolation 
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @param precip Logical indicating if o, p, s is precipitation data.
#' @param pr.threshold Integer. The minimum value that is considered as a non-zero precipitation. 
#' @param detrend logical. Detrend data prior to bias correction? Default. TRUE.
#' @param n.quantiles  Integer. Maximum number of quantiles to estimate. Default: same as data length.
#' @details DQM method developed by A. Canon, from \url{https://github.com/pacificclimate/ClimDown}, \url{https://cran.r-project.org/web/packages/ClimDown/}.
#'
#'  See also Cannon, A.J., S.R. Sobie, and T.Q. Murdock (2015) Bias Correction of GCM Precipitation by Quantile Mapping: How Well Do Methods Preserve Changes in Quantiles and Extremes?. J. Climate, 28, 6938–6959, \url{https://doi.org/10.1175/JCLI-D-14-00754.1}
#' @keywords internal
#' @author A. Cannon (acannon@uvic.ca), A. Casanueva

dqm <- function(o, p, s, precip, pr.threshold, n.quantiles, detrend=TRUE){
  
    if (all(is.na(o)) | all(is.na(p)) | all(is.na(s))) {
      return(yout=rep(NA, length(s)))
    
    } else{
      
      if(precip){
        # For ratio data, treat exact zeros as left censored values less than pr.threshold
        epsilon <- .Machine$double.eps
        o[o < pr.threshold & !is.na(o)] <- runif(sum(o < pr.threshold, na.rm=TRUE), min=epsilon, max=pr.threshold)
        p[p < pr.threshold & !is.na(p)] <- runif(sum(p < pr.threshold, na.rm=TRUE), min=epsilon, max=pr.threshold)
        s[s < pr.threshold & !is.na(s)] <- runif(sum(s < pr.threshold, na.rm=TRUE), min=epsilon, max=pr.threshold)
      }
      
      o.mn <- mean(o, na.rm=T)
      p.mn <- mean(p, na.rm=T)
      if(precip){
        s <- s/p.mn*o.mn
      } else{
        s <- s-p.mn+o.mn
      }
      
      if(detrend){
        s.mn <- lm.fit(cbind(1, seq_along(s)), s)$fitted
      } else{
        s.mn <- o.mn
      }
      if(is.null(n.quantiles)) n.quantiles <- max(length(o), length(p))
      tau <- c(0, (1:n.quantiles)/(n.quantiles+1), 1)
      if(precip & any(o < sqrt(.Machine$double.eps), na.rm=TRUE)){
        x <- quantile(p/p.mn, tau, na.rm=T)
        y <- quantile(o/o.mn, tau, na.rm=T)
        yout <- approx(x, y, xout=s/s.mn, rule=2:1)$y # if rule = 1, NAs are returned outside the training interval; if rule= 2, the value at the closest data extreme is used. rule = 2:1, if the left and right side extrapolation should differ.
        extrap <- is.na(yout)
        yout[extrap] <- max(o/o.mn, na.rm=T)*((s/s.mn)[extrap]/max(p/p.mn, na.rm=T)) # extrapolation on the upper tail
        yout <- yout*s.mn
        #yout.h <- approx(x, y, xout=p/p.mn, rule=1)$y*o.mn
      } else if(precip & !any(o < sqrt(.Machine$double.eps), na.rm=TRUE)){
        x <- quantile(p/p.mn, tau, na.rm=T)
        y <- quantile(o/o.mn, tau, na.rm=T)
        yout <- approx(x, y, xout=s/s.mn, rule=1)$y
        extrap.lower <- is.na(yout) & ((s/s.mn) < min(p/p.mn, na.rm=T))
        extrap.upper <- is.na(yout) & ((s/s.mn) > max(p/p.mn, na.rm=T))
        yout[extrap.lower] <- min(o/o.mn, na.rm=T)*((s/s.mn)[extrap.lower]/
                                                 min(p/p.mn, na.rm=T))
        yout[extrap.upper] <- max(o/o.mn, na.rm=T)*((s/s.mn)[extrap.upper]/
                                                 max(p/p.mn, na.rm=T))
        yout <- yout*s.mn
        #yout.h <- approx(x, y, xout=p/p.mn, rule=1)$y*o.mn
      } else{
        x <- quantile(p-p.mn, tau, na.rm=T)
        y <- quantile(o-o.mn, tau, na.rm=T)
        yout <- approx(x, y, xout=s-s.mn, rule=1)$y
        extrap.lower <- is.na(yout) & ((s-s.mn) < min(p-p.mn, na.rm=T))
        extrap.upper <- is.na(yout) & ((s-s.mn) > max(p-p.mn, na.rm=T))
        yout[extrap.lower] <- min(o-o.mn) + ((s-s.mn)[extrap.lower]-
                                                   min(p-p.mn, na.rm=T))
        yout[extrap.upper] <- max(o-o.mn) + ((s-s.mn)[extrap.upper]-
                                                   max(p-p.mn, na.rm=T))
        yout <- yout+s.mn
        #yout.h <- approx(x, y, xout=p-p.mn, rule=1)$y+o.mn
      }
      if(precip){
        yout[which(yout < sqrt(.Machine$double.eps))] <- 0
        #yout.h[which(yout.h < sqrt(.Machine$double.eps))] <- 0
      }
    return(yout)
    }
}


#' @title Quantile delta mapping
#' @description Quantile delta mapping 
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @param precip Logical indicating if o, p, s is precipitation data.
#' @param pr.threshold Integer. The minimum value that is considered as a non-zero precipitation. 
#' \code{precip = FALSE}.
#' @param jitter.factor Integer. Jittering to accomodate ties. Default: 0.01.
#' @param n.quantiles  Integer. Maximum number of quantiles to estimate. Default: same as data length.
#' @details QDM method developed by A. Canon, from \url{https://github.com/pacificclimate/ClimDown}, \url{https://cran.r-project.org/web/packages/ClimDown/}.
#' 
#' See also Cannon, A.J., S.R. Sobie, and T.Q. Murdock (2015) Bias Correction of GCM Precipitation by Quantile Mapping: How Well Do Methods Preserve Changes in Quantiles and Extremes?. J. Climate, 28, 6938–6959, \url{https://doi.org/10.1175/JCLI-D-14-00754.1}
#' @keywords internal
#' @author A. Cannon (acannon@uvic.ca), A. Casanueva

qdm <- function(o, p, s, precip, pr.threshold, n.quantiles, jitter.factor=0.01){
 
  # tau.s = F.s(x.s)
  # delta = x.s {/,-} F.p^-1(tau.s)
  # yout = F.o^-1(tau.s) {*,+} delta
  
  if (all(is.na(o)) | all(is.na(p)) | all(is.na(s))) {
    return(yout=rep(NA, length(s)))
  } else{
  
    # Apply a small amount of jitter to accomodate ties due to limited measurement precision
    o <- jitter(o, jitter.factor)
    p <- jitter(p, jitter.factor)
    s <- jitter(s, jitter.factor)
    
    # For ratio data, treat exact zeros as left censored values less than pr.threshold
    if(precip){
      epsilon <- .Machine$double.eps
      o[o < pr.threshold & !is.na(o)] <- runif(sum(o < pr.threshold, na.rm=TRUE), min=epsilon, max=pr.threshold)
      p[p < pr.threshold & !is.na(p)] <- runif(sum(p < pr.threshold, na.rm=TRUE), min=epsilon, max=pr.threshold)
      s[s < pr.threshold & !is.na(s)] <- runif(sum(s < pr.threshold, na.rm=TRUE), min=epsilon, max=pr.threshold)
    }
    
    # Calculate empirical quantiles using Weibull plotting position
    n <- max(length(o), length(p), length(s))
    if(is.null(n.quantiles)) n.quantiles <- n
    tau <- seq(1/(n+1), n/(n+1), length=n.quantiles)
    quant.o <- quantile(o, tau, type=6, na.rm=TRUE)
    quant.p <- quantile(p, tau, type=6, na.rm=TRUE)
    quant.s <- quantile(s, tau, type=6, na.rm=TRUE)
    
    # Apply QDM bias correction
    tau.s <- approx(quant.s, tau, s, rule=2)$y    
    if(precip){
      delta <- s/approx(tau, quant.p, tau.s, rule=2)$y # if rule= 2, the value at the closest data extreme is used
      yout <- approx(tau, quant.o, tau.s, rule=2)$y*delta
    } else{
      delta <- s - approx(tau, quant.p, tau.s, rule=2)$y
      yout <- approx(tau, quant.o, tau.s, rule=2)$y + delta
    }
    #yout.h <- approx(quant.p, quant.o, p, rule=2)$y
    
    # For precip data, set values less than threshold to zero
    if(precip){
      #yout.h[yout.h < pr.threshold] <- 0
      yout[yout < pr.threshold] <- 0
    }
    
    return(yout)
  }
}



#     biasCorrection.chunk.R Bias correction methods
#
#     Copyright (C) 2017 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Bias correction methods applied to each latitude (chunk)
#' @description Internal function that applies function \code{biasCorrection} 
#' to each latitude in \code{y} (chunk). 
#'
#' @template templateObsPredSim
#' @param output.path Path to the directory where the bias corrected *.rda objects are stored.
#' Default is the working directory.
#' @param method method applied. Current accepted values are \code{"eqm"}, \code{"delta"},
#'  \code{"scaling"}, \code{"pqm"} and \code{"gpqm"} \code{"variance"},\code{"loci"} and \code{"ptr"}. See details.
#' @param precipitation Logical for precipitation data (default to FALSE). If TRUE Adjusts precipitation 
#' frequency in 'x' (prediction) to the observed frequency in 'y'. This is a preprocess to bias correct 
#' precipitation data following Themeßl et al. (2012). To adjust the frequency, 
#' parameter \code{wet.threshold} is used (see below).
#' @param cross.val Logical (default to FALSE). Should cross-validation be performed? methods available are 
#' leave-one-out ("loo") and k-fold ("kfold") on an annual basis. The default option ("none") does not 
#' perform cross-validation.
#' @param folds Only requiered if \code{cross.val = "kfold"}. A list of vectors, each containing the years 
#' to be grouped in the corresponding fold.
#' @param wet.threshold The minimum value that is considered as a non-zero precipitation. Ignored when 
#' \code{precipitation = FALSE}. Default to 1 (assuming mm). See details on bias correction for precipitation.
#' @param window vector of length = 2 (or 1) specifying the time window width used to calibrate and the 
#' target days (days that are being corrected). The window is centered on the target day/s 
#' (window width >= target days). Default to \code{NULL}, which considers the whole period.
#' @param scaling.type Character indicating the type of the scaling method. Options are \code{"additive"} (default)
#' or \code{"multiplicative"} (see details). This argument is ignored if \code{"scaling"} is not 
#' selected as the bias correction method.
#' @param  fitdistr.args Further arguments passed to function \code{\link[MASS]{fitdistr}} 
#' (\code{densfun}, \code{start}, \code{...}). Only used when applying the "pqm" method 
#' (parametric quantile mapping). Please, read the \code{\link[MASS]{fitdistr}} help 
#' document  carefully before setting the parameters in \code{fitdistr.args}.
#' @param n.quantiles Integer indicating the number of quantiles to be considered when method = "eqm". Default is NULL, 
#' that considers all quantiles, i.e. \code{n.quantiles = length(x[i,j])} (being \code{i} and \code{j} the 
#' coordinates in a single location).
#' @param extrapolation Character indicating the extrapolation method to be applied to correct values in  
#' \code{newdata} that are out of the range of \code{x}. Extrapolation is applied only to the \code{"eqm"} method, 
#' thus, this argument is ignored if other bias correction method is selected. Default is \code{"none"} (do not extrapolate).
#' @param theta numeric indicating  upper threshold (and lower for the left tail of the distributions, if needed) 
#' above which precipitation (temperature) values are fitted to a Generalized Pareto Distribution (GPD). 
#' Values below this threshold are fitted to a gamma (normal) distribution. By default, 'theta' is the 95th 
#' percentile (5th percentile for the left tail). Only for \code{"gpqm"} method.
#' @param join.members Logical indicating whether members should be corrected independently (\code{FALSE}, the default),
#'  or joined before performing the correction (\code{TRUE}). It applies to multimember grids only (otherwise ignored).
#' @template templateParallelParams
#'  
#' @details
#' 
#' The methods available are \code{"eqm"}, \code{"delta"}, 
#' \code{"scaling"}, \code{"pqm"}, \code{"gpqm"}\code{"loci"}, 
#' \code{"ptr"}  (the four latter used only for precipitation) and 
#' \code{"variance"} (only for temperature).
#' 
#'  These are next briefly described: 
#'  
#' \strong{Delta}
#'
#' This method consists on adding to the observations the mean change signal (delta method).
#' This method is applicable to any kind of variable but it is preferable to avoid it for bounded variables
#' (e.g. precipitation, wind speed, etc.) because values out of the variable range could be obtained
#' (e.g. negative wind speeds...). This method corresponds to case g=1 and f=0 in Amengual et al. 2012. 
#' 
#' \strong{Scaling}
#' 
#' This method consists on scaling the simulation  with the difference (additive) or quotient (multiplicative) 
#' between the observed and simulated means in the train period. The \code{additive} or \code{multiplicative}
#' correction is defined by parameter \code{scaling.type} (default is \code{additive}).
#' The additive version is preferably applicable to unbounded variables (e.g. temperature) 
#' and the multiplicative to variables with a lower bound (e.g. precipitation, because it also preserves the frequency). 
#' 
#' 
#' \strong{eqm}
#' 
#' Empirical Quantile Mapping. This is a very extended bias correction method which consists on calibrating the simulated Cumulative Distribution Function (CDF) 
#' by adding to the observed quantiles both the mean delta change and the individual delta changes in the corresponding quantiles. 
#' This is equivalent to f=g=1 in Amengual et al. 2012. This method is applicable to any kind of variable.
#' 
#' 
#' \strong{pqm}
#' 
#' Parametric Quantile Mapping. It is based on the initial assumption that both observed and simulated intensity distributions are well approximated by a given distribution
#' (see \code{\link[MASS]{fitdistr}} to check available distributions), therefore is a parametric q-q map that uses the theorical instead of the empirical distribution.
#' For instance, the gamma distribution is described in Piani et al. 2010 and is applicable to precipitation. Other example is the weibull distribution, which
#' is applicable to correct wind data (Tie et al. 2014).
#'  
#' \strong{gpqm}
#'  
#' Generalized Quantile Mapping (described in Gutjahr and Heinemann 2013) is also a parametric quantile mapping (see
#' method 'pqm') but using two teorethical distributions, the gamma distribution and Generalized Pareto Distribution (GPD).
#' By default, It applies a gamma distribution to values under the threshold given by the 95th percentile 
#' (following Yang et al. 2010) and a general Pareto distribution (GPD) to values above the threshold. the threshold above 
#' which the GPD is fitted is the 95th percentile of the observed and the predicted wet-day distribution, respectively. 
#' The user can specify a different threshold by modifying the parameter theta. It is applicable to precipitation data. 
#' 
#' 
#' \strong{variance}
#' 
#' Variance scaling of temperature. This method is described in Chen et al. 2011. It is applicable only to temperature. It corrects
#' the mean and variance of temperature time series.
#' 
#' \strong{loci}
#' 
#' Local intensity scaling of precipitation. This method is described in Schmidli et al. 2006. It adjust the mean as well as both wet-day frequencies and wet-day intensities.
#' The precipitation threshold is calculated such that the number of simulated days exceeding this threshold matches the number of observed days with precipitation larger than 1 mm.
#' 
#'\strong{ptr}
#'
#' Power transformation of precipitation. This method is described in Leander and Buishand 2007 and is applicable only to precipitation. It adjusts the variance statistics of precipitation
#' time series in an exponential form. The power parameter is estimated on a monthly basis using a 90-day window centered on the interval. The power is defined by matching the coefficient
#' of variation of corrected daily simulated precipitation with the coefficient of variation of observed daily precipitation. It is calculated by root-finding algorithm using Brent's method.
#'
#'
#' @section Note on the bias correction of precipitation:
#' 
#' In the case of precipitation a frequency adaptation has been implemented in all versions of 
#' quantile mapping to alleviate the problems arising when the dry day frequency in the raw model output is larger
#'  than in the observations (Wilcke et al. 2013). 
#'  
#'  The precipitation subroutines are switched-on when the variable name of the grid 
#'  (i.e., the value returned by \code{gridData$Variable$varName}) is one of the following: 
#'  \code{"pr"}, \code{"tp"} (this is the standard name defined in the vocabulary (\code{\link[loadeR]{C4R.vocabulary}}), \code{"precipitation"} or \code{"precip"}.
#'  Thus, caution must be taken to ensure that the correct bias correction is being undertaken when dealing with
#'  non-standard variables.
#'     
#' 
#' @seealso \code{\link{isimip}} for a trend-preserving method of model calibration and \code{\link{quickDiagnostics}} 
#' for an outlook of the results.
#' @return Calibrated grids for each latitudinal chunk (with the same spatio-temporal extent than the chunked input \code{"y"}). 
#' This objects are saved in the specified \code{output.path}. The object obatained in the workspace 
#' is a charecter string of the listed files.
#' @family downscaling
#' 
#' @importFrom transformeR redim subsetGrid getYearsAsINDEX getDim
#' @importFrom abind adrop
#'
#' @references
#'
#' \itemize{
#' \item R.A.I. Wilcke, T. Mendlik and A. Gobiet (2013) Multi-variable error correction of regional climate models. Climatic Change, 120, 871-887
#'
#' \item A. Amengual, V. Homar, R. Romero, S. Alonso, and C. Ramis (2012) A Statistical Adjustment of Regional Climate Model Outputs to Local Scales: Application to Platja de Palma, Spain. J. Clim., 25, 939-957
#'
#' \item C. Piani, J. O. Haerter and E. Coppola (2009) Statistical bias correction for daily precipitation in regional climate models over Europe, Theoretical and Applied Climatology, 99, 187-192
#'
#' \item O. Gutjahr and G. Heinemann (2013) Comparing precipitation bias correction methods for high-resolution regional climate simulations using COSMO-CLM, Theoretical and Applied Climatology, 114, 511-529
#' 
#' \item M. R. Tye, D. B. Stephenson, G. J. Holland and R. W. Katz (2014) A Weibull Approach for Improving Climate Model Projections of Tropical Cyclone Wind-Speed Distributions, Journal of Climate, 27, 6119-6133
#' 
#' }
#' @author M. Iturbide
#' @examples {
#' # empirical
#' eqm1 <- biasCorrection.chunk(output.path = getwd(),
#'                              n.chunks = 10,
#'                              login.UDG,
#'                              dataset.y = "WATCH_WFDEI",
#'                              dataset.x = "CORDEX-EUR11_EC-EARTH_r12i1p1_historical_RCA4_v1",
#'                              dataset.newdata = "CORDEX-EUR11_EC-EARTH_r12i1p1_historical_RCA4_v1",
#'                              loadGridData.args = list(var = "tasmax", 
#'                                                       years = 1971:2000, 
#'                                                       lonLim = c(-10, 10), 
#'                                                       latLim = c(36, 45)),
#'                              biasCorrection.args = list(precipitation = FALSE,
#'                                                         method = c("delta", "scaling", "eqm", "pqm", "gpqm", "loci"),
#'                                                         cross.val = c("none", "loo", "kfold"),
#'                                                         folds = NULL,
#'                                                         window = NULL,
#'                                                         scaling.type = c("additive", "multiplicative"),
#'                                                         fitdistr.args = list(densfun = "normal"),
#'                                                         wet.threshold = 1,
#'                                                         n.quantiles = NULL,
#'                                                         extrapolation = c("none", "constant"), 
#'                                                         theta = .95,
#'                                                         join.members = FALSE,
#'                                                         parallel = FALSE,
#'                                                         max.ncores = 16,
#'                                                         ncores = NULL))
#' }

# eqm1 <- biasCorrection.chunk(output.path = getwd(),
#                              n.chunks = 10,
#                              loginUDG.args = list(username = "miturbide", password = "iturbide.14"),
#                              dataset.y = "WATCH_WFDEI",
#                              dataset.x = "CORDEX-EUR11_EC-EARTH_r12i1p1_historical_RCA4_v1",
#                              dataset.newdata = "CORDEX-EUR11_EC-EARTH_r12i1p1_historical_RCA4_v1",
#                              loadGridData.args = list(var = "tasmax",
#                                                       years = 1971:2000,
#                                                       lonLim = c(-10, 10),
#                                                       latLim = c(36, 45)),
#                              biasCorrection.args = list(precipitation = FALSE,
#                                                         method = c("delta", "scaling", "eqm", "pqm", "gpqm", "loci"),
#                                                         cross.val = c("none", "loo", "kfold"),
#                                                         folds = NULL,
#                                                         window = NULL,
#                                                         scaling.type = c("additive", "multiplicative"),
#                                                         fitdistr.args = list(densfun = "normal"),
#                                                         wet.threshold = 1,
#                                                         n.quantiles = NULL,
#                                                         extrapolation = c("none", "constant"),
#                                                         theta = .95,
#                                                         join.members = FALSE,
#                                                         parallel = FALSE,
#                                                         max.ncores = 16,
#                                                         ncores = NULL))
# 
# biasCorrection.chunk <- function(output.path = getwd(),
#                            n.chunks = 10,
#                            loginUDG.args = list(NULL),
#                            dataset.y = "WATCH_WFDEI",
#                            dataset.x = "CORDEX-EUR11_EC-EARTH_r12i1p1_historical_RCA4_v1",
#                            dataset.newdata = "CORDEX-EUR11_EC-EARTH_r12i1p1_historical_RCA4_v1",
#                            loadGridData.args = list(var = "tasmax", 
#                                                     years = 1971:2000, 
#                                                     lonLim = c(-10, 10), 
#                                                     latLim = c(36, 45)),
#                            biasCorrection.args = list(precipitation = FALSE,
#                               method = c("delta", "scaling", "eqm", "pqm", "gpqm", "loci"),
#                               cross.val = c("none", "loo", "kfold"),
#                               folds = NULL,
#                               window = NULL,
#                               scaling.type = c("additive", "multiplicative"),
#                               fitdistr.args = list(densfun = "normal"),
#                               wet.threshold = 1,
#                               n.quantiles = NULL,
#                               extrapolation = c("none", "constant"), 
#                               theta = .95,
#                               join.members = FALSE,
#                               parallel = FALSE,
#                               max.ncores = 16,
#                               ncores = NULL)) {
#       suppressWarnings(dir.create(output.path))
#       message("[", Sys.time(), "] Rdata will be saved in ", output.path)
#       do.call("loginUDG", login.UDG)
#       di.y <- dataInventory(dataset.y)
#       lats.y <- di.y[[loadGridData.args[["var"]]]]$Dimensions$lat$Values
#       if (is.null(lats.y)) stop("dataset.y does not contain the requested variable")
#       lats.y <- lats.y[which.min(abs(lats.y - latLim[1]))[1]:(which.min(abs(lats.y - latLim[2]))[1] + 1)]
#       n.lats.y <- length(lats.y)
#       n.lat.chunk <- ceiling(n.lats.y/n.chunks)
#       aux.ind <- rep(1:(n.chunks - 1), each = n.lat.chunk)
#       ind <- c(aux.ind, rep((max(aux.ind) + 1), each = n.lats.y - length(aux.ind)))
#       message("[", Sys.time(), "] y contains ", n.lats.y, " latitudes. Bias correction will be applied in ",n.chunks, " chunks of about ", n.lat.chunk, " latitudes.")
#       lat.list <- split(lats.y, f = ind)
#       lat.range.chunk <- lapply(lat.list, range)
#       lat.range.chunk.x <- lapply(lat.range.chunk, function(x) c(x[1] - 3, x[2] + 3))
#       
#       file.dir <- character()
#       for (i in 1:length(lat.range.chunk)) {
#             loadGridData.args[["dataset"]] <- dataset.y
#             loadGridData.args[["latLim"]] <- lat.range.chunk[[i]]
#             y <- do.call("loadGridData", loadGridData.args)
#             loadGridData.args[["dataset"]] <- dataset.x
#             loadGridData.args[["latLim"]] <- lat.range.chunk.x[[i]]
#             x <- do.call("loadGridData", loadGridData.args)
#             if (!is.null(dataset.newdata)) {
#                   loadGridData.args[["dataset"]] <- dataset.newdata
#                   newdata <- do.call("loadGridData", loadGridData.args)
#             } else {
#                   newdata <- NULL
#             }
#             bc <- do.call("biasCorrection", biasCorrection.args)
#             file.dir[i] <- paste0(output.path, "/chunk", i, ".rda")
#             save(bc, file = paste0(output.path, "/chunk", i, ".rda"))
#       }
#       return(file.dir)
# }
# 
