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
#'  \code{"scaling"}, \code{"gqm"} and \code{"gpqm"} \code{"variance"},\code{"loci"} and \code{"ptr"}. See details.
#' @param precipitation Logical indicating if the data to be corrected is precipitation data.
#' @param cross.val Should cross-validation be performed? methods available are leave-one-out ("loo") 
#' and k-fold ("kfold") on an annual basis. The default option ("none") does not perform cross-validation.
#' @param folds Only requiered if cross.val = "kfold". A list of vectors, each containing the years to be grouped in 
#' the corresponding fold.
#' @param wet.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm). See details on bias correction for precipitation.
#' @param window vector of length = 2 specifying the time window width used to calibrate and the target days (days that are being corrected).
#'  The window is centered on the target day/s (window width >= target days). 
#' Default to \code{NULL}, which considers the whole period available.
#' 
#' @param scaling.type Character indicating the type of the scaling method. Options are \code{"additive"} (default)
#' or \code{"multiplicative"} (see details). This argument is ignored if \code{"scaling"} is not 
#' selected as the bias correction method.
#' @param n.quantiles Integer indicating the number of quantiles to be considered when method = "eqm". Default is NULL, 
#' that considers all quantiles, i.e. \code{n.quantiles = length(x[i,j])} (being \code{i} and \code{j} the coordinates in a single location).
#' @param extrapolation Character indicating the extrapolation method to be applied to correct values in  
#' \code{newdata} that are out of the range of \code{x}. Extrapolation is applied only to the \code{"eqm"} method, 
#' thus, this argument is ignored if other bias correction method is selected. Default is \code{"none"} (do not extrapolate).
#' @param theta numeric indicating  upper threshold (and lower for the left tail of the distributions, if needed) 
#' above which precipitation (temperature) values are fitted to a Generalized Pareto Distribution (GPD). 
#' Values below this threshold are fitted to a gamma (normal) distribution. By default, 'theta' is the 95th 
#' percentile (5th percentile for the left tail). Only for \code{"gpqm"} method.
#' @param join.members Logical indicating whether members should be corrected independently (\code{FALSE}, the default),
#'  or joined before performing the correction (\code{TRUE}). It applies to multimember grids only (otherwise ignored).
#'  
#' @details
#' 
#' The methods available are \code{"eqm"}, \code{"delta"}, 
#' \code{"scaling"}, \code{"gqm"}, \code{"gpqm"}\code{"loci"}, 
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
#' (e.g. negative wind speeds...).
#' 
#' \strong{Scaling}
#' 
#' This method consists on scaling the simulation  with the difference (additive) or quotient (multiplicative) 
#' between the observed and simulated means in the train period. The \code{additive} or \code{multiplicative}
#' correction is defined by parameter \code{scaling.type} (default is \code{additive}).
#' The additive version is preferably applicable to unbounded variables (e.g. temperature) 
#' and the multiplicative to variables with a lower bound (e.g. precipitation, because it also preserves the frequency). 
#' 
#' \strong{eqm}
#' 
#' Empirical Quantile Mapping. This is a very extended bias correction method which consists on calibrating the simulated Cumulative Distribution Function (CDF) 
#' by adding to the observed quantiles both the mean delta change and the individual delta changes in the corresponding quantiles. 
#' This method is applicable to any kind of variable.
#' 
#' \strong{gqm}
#' 
#' Gamma Quantile Mapping. This method is described in Piani et al. 2010 and is applicable only to precipitation. It is based on the initial assumption that both observed
#' and simulated intensity distributions are well approximated by the gamma distribution, therefore is a parametric q-q map 
#' that uses the theorical instead of the empirical distribution. 
#'  
#' \strong{gpqm}
#'  
#' Generalized Quantile Mapping. This method is described in Gutjahr and Heinemann 2013. It is applicable only to precipitation and is similar to the Piani method. It applies a 
#' gamma distribution to values under the threshold given by the 95th percentile (following Yang et al. 2010) and a general Pareto 
#' distribution (GPD) to values above the threshold.
#' 
#' 
#' \strong{variance}
#' 
#' Variance scaling of temperature. This method is described in Chen et al. 2011. It is applicable only to temperature. It corrects
#' the mean and variance of temperature time series.
#' 
#' \strong{loci}
#' 
#' Local intensity scaling of precipitation. This methode is described in Schmidli et al. 2006. It adjust the mean as well as both wet-day frequencies and wet-day intensities.
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
#' qqmap to alleviate the problems arising when the dry day frequency in the raw model output is larger
#'  than in the observations (Wilcke et al. 2013). 
#'  
#'  The precipitation subroutines are switched-on when the variable name of the grid 
#'  (i.e., the value returned by \code{gridData$Variable$varName}) is one of the following: 
#'  \code{"pr"}, \code{"tp"} (this is the standard name defined in the vocabulary (\code{\link[loadeR]{UDG.vocabulary}}), \code{"precipitation"} or \code{"precip"}.
#'  Thus, caution must be taken to ensure that the correct bias correction is being undertaken when dealing with
#'  non-standard variables.
#'     
#' 
#' @seealso \code{\link{isimip}} for a trend-preserving method of model calibration and \code{\link{quickDiagnostics}} 
#' for an outlook of the results.
#' @return A calibrated grid of the same spatio-temporal extent than the input \code{"y"}
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
#' }
#' @author S. Herrera, M. Iturbide, J. Bedia
#' @export
#' @examples {
#' data("VALUE_Iberia_pr")
#' data("NCEP_Iberia_pr")
#' y <- VALUE_Iberia_pr
#' x <- NCEP_Iberia_pr
#' x$Data <- x$Data*86400
#' 
#' eqm1 <- biasCorrection(y = y, x = x, 
#'                        method = "eqm",
#'                        window = NULL)
#' eqm1win <- biasCorrection(y = y, x = x, 
#'                           method = "eqm",
#'                           extrapolation = "none",
#'                           window = c(31, 1))
#' eqm1folds <- biasCorrection(y = y, x = x,
#'                             method = "eqm",
#'                             window = c(31, 1),
#'                             cross.val = "kfold",
#'                             folds = list(1983:1989, 1990:1996, 1997:2002))
#' 
#' quickDiagnostics(y, x, eqm1)
#' quickDiagnostics(y, x, eqm1, location = c(-2, 43))
#' quickDiagnostics(y, x, eqm1win, location = c(-2, 43))
#' quickDiagnostics(y, x, eqm1folds, location = c(-2, 43))
#' }



biasCorrection <- function(y, x, newdata = NULL, precipitation = FALSE,
                           method = c("delta", "scaling", "eqm", "gqm", "gpqm", "loci"),
                           cross.val = c("none", "loo", "kfold"),
                           folds = NULL,
                           window = NULL,
                           scaling.type = c("additive", "multiplicative"),
                           wet.threshold = 1,
                           n.quantiles = NULL,
                           extrapolation = c("none", "constant"), 
                           theta = .95,
                           join.members = FALSE) {
      method <- match.arg(method, choices = c("delta", "scaling", "eqm", "gqm", "gpqm", "loci", "ptr", "variance"))
      cross.val <- match.arg(cross.val, choices = c("none", "loo", "kfold"))
      scaling.type <- match.arg(scaling.type, choices = c("additive", "multiplicative"))
      extrapolation <- match.arg(extrapolation, choices = c("none", "constant"))
      stopifnot(is.logical(join.members))
      nwdatamssg <- TRUE
      if (is.null(newdata)) {
            newdata <- x 
            nwdatamssg <- FALSE
      }
      if ("loc" %in% getDim(y) & isTRUE(join.members)) {
            join.members <- FALSE
            warning("The option 'join.members=TRUE' is currently not supported for station data predictand. It was reset to 'FALSE'.")
      }
      if (cross.val == "none") {
            output <- biasCorrectionXD(y = y, x = x, newdata = newdata, 
                                       precipitation = precipitation,
                                       method = method,
                                       window = window,
                                       scaling.type = scaling.type,
                                       pr.threshold = wet.threshold, 
                                       n.quantiles = n.quantiles, 
                                       extrapolation = extrapolation, 
                                       theta = theta,
                                       join.members = join.members)
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
                        newdata2 <- subsetGrid(x, years = target.year)
                        xx <- subsetGrid(x, years = rest.years)
                        message("Validation ", i, ", ", length(unique(years)) - i, " remaining")
                        biasCorrectionXD(y = yy, x = xx, newdata = newdata2, precipitation = precipitation,
                                          method = method,
                                          window = window,
                                          scaling.type = scaling.type,
                                          pr.threshold = wet.threshold, n.quantiles = n.quantiles, extrapolation = extrapolation, 
                                          theta = theta, join.members = join.members)
                  })
                  al <- which(getDim(x) == "time")
                  Data <- sapply(output.list, function(n) unname(n$Data), simplify = FALSE)
                  bindata <- unname(do.call("abind", c(Data, along = al)))
                  output <- output.list[[1]]
                  dimNames <- attr(output$Data, "dimensions")
                  output$Data <- bindata
                  attr(output$Data, "dimensions") <- dimNames
                  output$Dates <- x$Dates
      }
      return(output)
}

#' @keywords internal
#' @importFrom transformeR redim subsetGrid getDim

biasCorrectionXD <- function(y, x, newdata, precipitation, 
                             method = method,
                             window = NULL,
                             scaling.type = c("additive", "multiplicative"),
                             pr.threshold = 1, n.quantiles = NULL, extrapolation = c("none", "constant"), 
                             theta = .95,
                             join.members = join.members) {
      obso <- y
      pred <- x
      sim <- newdata
      delta.method <- method == "delta"
      precip <- precipitation
      message("[", Sys.time(), "] Argument precipitation is set as ", precip, ", please ensure that this matches your data.")
      if ("loc" %in% getDim(obso)) {
            station <- TRUE
            obs <- redim(obso, member = FALSE)
            x <- obs$xyCoords[,1]
            y <- obs$xyCoords[,2]
            ito <- which(getDim(obs) == "time")
            ind <- cbind(1:dim(obs$Data)[2], rep(1, dim(obs$Data)[2]), 1:dim(obs$Data)[2])
      } else {
            station <- FALSE
            x <- obso$xyCoords$x
            y <- obso$xyCoords$y
            obs <- obso
            ind1 <- expand.grid(1:length(y), 1:length(x))
            ind <- cbind(ind1, ind1[,2])
      }
      bc <- obso
      if (isTRUE(join.members)) {
            pred <- flatMemberDim(pred)
            sim <- flatMemberDim(sim)
      }
      pred <- redim(pred, member = TRUE, runtime = TRUE)
      sim <- if (!is.null(sim)) {
            redim(sim, member = TRUE, runtime = TRUE)
      } else {
            pred
      }
      ito <- which(getDim(obs) == "time")
      n.run <- dim(sim$Data)[1]
      n.mem <- dim(sim$Data)[2]
      if (delta.method) {
            run <- array(dim = c(1, n.mem, dim(obs$Data)))
      } else {
            its <- which(getDim(sim) == "time")
            run <- array(dim = c(1, n.mem, dim(sim$Data)[its], dim(obs$Data)[-ito]))
      }
      lrun <- lapply(1:n.run, function(k) {    # loop for runtimes
            pre <- subsetGrid(pred, runtime = k, drop = FALSE)
            si <- subsetGrid(sim, runtime = k, drop = FALSE)
            if (delta.method) {
                  mem <- array(dim = c(1, 1, dim(obs$Data)))
            } else {
                  mem <- array(dim = c(1, 1, dim(sim$Data)[its], dim(obs$Data)[-ito]))
            }
            lmem <- lapply(1:n.mem, function(l) { # loop for members
                  p <- subsetGrid(pred, members = l, drop = FALSE)
                  s <- subsetGrid(sim, members = l, drop = FALSE)
                  #join members
                  if (!isTRUE(join.members)) {
                        message("[", Sys.time(), "] Bias-correcting member ", l, " out of ", n.mem, "...")
                  } else {
                        message("[", Sys.time(), "] Bias-correcting ", attr(pred, "orig.mem.shape"), " members considering their joint distribution...")
                  }
                  #window
                  if (!is.null(window)) {
                        win <- getWindowIndex(y = obs, newdata = s, window = window, delta.method = delta.method)
                  } else {
                        win <- list()
                        win[["Window1"]] <- list("window" = 1:getShape(obs)["time"], "step" = 1:getShape(s)["time"])
                  }
                  message("[", Sys.time(), "] Number of windows considered: ", length(win), "...")
                  #correcting locations
                  for (i in 1:nrow(ind)) {
                        # Apply bias correction methods
                        for (j in 1:length(win)) {
                              yind <- win[[j]]$window
                              outind <- win[[j]]$step
                              if (delta.method) {
                                    yind <- win[[j]]$deltaind
                                    outind <- win[[j]]$deltaind
                              } 
                              suppressWarnings(
                              mem[,,outind,ind[i,1],ind[i,2]] <- biasCorrection1D(obs$Data[yind,ind[i,1],ind[i,2]],
                                                                            subsetGrid(p, latLim = y[ind[i,1]], lonLim = x[ind[i,3]], outside = TRUE, drop = FALSE)$Data[1,1,win[[j]]$window,1,1],
                                                                            subsetGrid(s, latLim = y[ind[i,1]], lonLim = x[ind[i,3]], outside = TRUE, drop = FALSE)$Data[1,1,win[[j]]$step,1,1], 
                                                                            method = method,
                                                                            scaling.type = scaling.type,
                                                                            precip = precip,
                                                                            pr.threshold = pr.threshold,
                                                                            n.quantiles = n.quantiles,
                                                                            extrapolation = extrapolation,
                                                                            theta = theta)
                              )
                              
                        }
                  }
                  return(mem)
            }) 
            run[,,,,] <- abind(lmem, along = 2)
            return(run)
      }) #end loop for runtimes
      bc$Data <- unname(abind(lrun, along = 1))
      attr(bc$Data, "dimensions") <- attr(sim$Data, "dimensions")
      ## Recover the member dimension when join.members=TRUE:
      if (isTRUE(join.members)) {
            bc <- recoverMemberDim(pred, bc, newdata)
      } else {
            bc$Dates <- sim$Dates
            bc$InitializationDates <- sim$InitializationDates
            bc$Members <- sim$Members
      }
      attr(bc$Variable, "correction") <- method
      bc <- redim(bc, drop = TRUE)
      if(station & !"loc" %in% bc) bc <- redim(bc, member = FALSE, loc = TRUE)
      message("[", Sys.time(), "] Done.")
      return(bc)
}


#' @title Get the index of the window days.
#' @description Get the index of the days that corresponding to the window and the target days centered in it.
#' 
#' @param y A grid or station data containing the observed climate data for the training period
#' @param newdata A grid containing the simulated climate for the test period.
#' @param method method applied. Current accepted values are \code{"eqm"}, \code{"delta"},
#'  \code{"scaling"}, \code{"gqm"} and \code{"gpqm"} \code{"variance"},\code{"loci"} and \code{"ptr"}. See details.
#' @param window vector of length = 2 specifying the time window width used to calibrate and the target days (days that are being corrected).
#'  The window is centered on the target day/s (window width >= target days). 
#' @param delta.method logical (default is FALSE)
#' @keywords internal
#' @author M. Iturbide

getWindowIndex <- function(y, newdata, window, delta.method = FALSE){
      step <- window[2]
      window <- window[1]
      if (window - step < 0) stop("The first argument of window must be equal or higher than the second. See ?biasCorrection")
      datesList <- as.POSIXct(y$Dates$start, tz = "GMT", format = "%Y-%m-%d")
      yearList <- unlist(strsplit(as.character(datesList), "[-]"))
      dayListObs <- array(data = c(as.numeric(yearList[seq(2,length(yearList),3)]),as.numeric(yearList[seq(3,length(yearList),3)])), dim = c(length(datesList),2))
      dayList <- unique(dayListObs,index.return = datesList)
      annual <- TRUE
      if (nrow(dayList) < 360) annual <- FALSE
      indDays <- array(data = NaN, dim = c(length(datesList),1))
      for (d in 1:dim(dayList)[1]) {
            indDays[which(sqrt((dayListObs[,1] - dayList[d,1]) ^ 2 + (dayListObs[,2] - dayList[d,2]) ^ 2) == 0)] <- d
      }
      datesList <- as.POSIXct(newdata$Dates$start, tz = "GMT", format = "%Y-%m-%d")
      yearList <- unlist(strsplit(as.character(datesList), "[-]"))
      dayListSim <- array(data = c(as.numeric(yearList[seq(2,length(yearList),3)]),as.numeric(yearList[seq(3,length(yearList),3)])), dim = c(length(datesList),2))
      indDaysSim <- array(data = NaN, dim = c(length(datesList),1))
      for (d in 1:dim(dayList)[1]) {
            indDaysSim[which(sqrt((dayListSim[,1] - dayList[d,1]) ^ 2 + (dayListSim[,2] - dayList[d,2]) ^ 2) == 0)] <- d
      }
      steps <- floor(dim(dayList)[1]/step)
      #steps loop
      output <- list()
      for (j in 1:steps) {
            days <- ((j - 1) * step + 1):((j - 1) * step + step)
            if (j == steps) days <- days[1]:dim(dayList)[1]
            indObs <- lapply(1:length(days), function(h){
                  which(indDays == days[h])
            })
            indObs <- sort(do.call("abind", indObs))
            head <- floor((window - step)/2)
            tail <- head
            before <- after <- FALSE
            if (!annual) {
                  before <- min(indDays[indObs]) - 1 - head < 1
                  if (before) head <- head + (min(indDays[indObs]) - 1 - head) 
                  after <- max(indDays[indObs]) + tail > nrow(dayList)
                  if (after) tail <- nrow(dayList) - max(indDays[indObs])   
            }
            indObsWindow <- array(data = NA, dim = c((head + step + tail)*length(indObs)/step,1))
            breaks <- c(which(diff(indObs) != 1), length(indObs))
            for (d in 1:length(breaks)) {
                  if (d == 1) {
                        piece <- indObs[1:breaks[1]]
                  } else {
                        piece <- indObs[(breaks[d - 1] + 1):breaks[d]]
                  }
                  suppressWarnings(indObsWindow[((d - 1) * (head + step + tail) + 1):(d * (head + step + tail))] <- 
                                         (min(piece, na.rm = TRUE) - head):(max(piece, na.rm = TRUE) + tail))
            }
            if (annual) {
                  indObsWindow[which(indObsWindow <= 0)] <- 1
                  indObsWindow[which(indObsWindow >  length(indDays))] <- length(indDays)
                  indObsWindow <- unique(indObsWindow)
            }
            indSim <- lapply(1:length(days), function(h){
                  which(indDaysSim == days[h])
            })
            indSim <- sort(do.call("abind", indSim))
            names(indSim) <- newdata$Dates$start[indSim]
            indObsWindow <- indObsWindow[which(!is.na(y$Dates$start[indObsWindow]))]
            names(indObsWindow) <- y$Dates$start[indObsWindow]
            names(indObs) <- y$Dates$start[indObs]
            output[[paste0("Window", j)]] <- list("window" = indObsWindow, "step" = indSim)
            if (delta.method) output[[paste0("Window", j)]][["deltaind"]] <- indObs
      }
      return(output)
}


#' @title Bias correction methods on 1D data
#' @description Implementation of several standard bias correction methods
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @param method method applied. Current accepted values are \code{"eqm"}, \code{"delta"},
#'  \code{"scaling"}, \code{"gqm"}, \code{"gpqm"}, \code{"variance"}, \code{"loci"} and \code{"ptr"}. 
#' @param scaling.type Character indicating the type of the scaling method. Options are \code{"additive"} (default)
#' or \code{"multiplicative"} (see details). This argument is ignored if \code{"scaling"} is not selected as the bias correction method.
#' @param precip Logical indicating if o, p, s is precipitation data.
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm). 
#' @param extrapolation Character indicating the extrapolation method to be applied to correct values in  
#' \code{"s"} that are out of the range of \code{"p"}. Extrapolation is applied only to the \code{"eqm"} method, 
#' thus, this argument is ignored if other bias correction method is selected. Default is \code{"none"} (do not extrapolate).
#' @param theta numeric indicating  upper threshold (and lower for the left tail of the distributions, if needed) 
#' above which precipitation (temperature) values are fitted to a Generalized Pareto Distribution (GPD). 
#' Values below this threshold are fitted to a gamma (normal) distribution. By default, 'theta' is the 95th 
#' percentile (5th percentile for the left tail). Only for \code{"gpqm"} method.
#' @keywords internal
#' @author M. Iturbide

biasCorrection1D <- function(o, p, s,
                             method = method, 
                             scaling.type = scaling.type,
                             precip = FALSE, 
                             pr.threshold = 1,
                             n.quantiles = NULL,
                             extrapolation = extrapolation, 
                             theta = .95) {
      if (method == "delta") {
            delta(o, p, s)
      } else if (method == "scaling") {
            scaling(o, p, s, scaling.type)
      } else if (method == "eqm") {
            eqm(o, p, s, precip, pr.threshold, n.quantiles, extrapolation)
      } else if (method == "gqm") {
            gqm(o, p, s, precip, pr.threshold)
      } else if (method == "gpqm") {
            gpqm(o, p, s, precip, pr.threshold, theta)
      } else if (method == "variance") {
            variance(o, p, s, precip)
      } else if (method == "loci") {
            loci(o, p, s, precip, pr.threshold)
      } else if (method == "ptr") {
            ptr(o, p, s, precip)
      }
}


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
            s - mean(p) + mean(o)
      } else if (scaling.type == "multiplicative") {
            (s/mean(p)) * mean(o)
      }
}



#' @title Gamma Quantile Mapping method for bias correction
#' @description Implementation of Gamma Quantile Mapping method for bias correction 
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{x}, but considering the test period.
#' @param precip Logical indicating if o, p, s is precipitation data.
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm). 
#' @importFrom MASS fitdistr
#' @importFrom stats pgamma qgamma
#' @keywords internal
#' @author S. Herrera and M. Iturbide

gqm <- function(o, p, s, precip, pr.threshold){
      if (precip == FALSE) {
            stop("method gqm is only applied to precipitation data")
      } else {
            threshold <- pr.threshold
            if (any(!is.na(o))) {
                  params <-  norain(o, p, threshold)
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
                  obsGamma <-  tryCatch({fitdistr(o[ind],"gamma")}, error = function(err){stop("There are not precipitation days in y for the window length selected in one or more locations. Try to enlarge the window")})
                  ind <- which(p > 0 & !is.na(p))
                  prdGamma <- tryCatch({fitdistr(p[ind],"gamma")}, error = function(err){stop("There are not precipitation days in x for the window length selected in one or more locations. Try to enlarge the window")})
                  rain <- which(s > Pth & !is.na(s))
                  noRain <- which(s <= Pth & !is.na(s))
                  auxF <- pgamma(s[rain], prdGamma$estimate[1], rate = prdGamma$estimate[2])
                  s[rain] <- qgamma(auxF, obsGamma$estimate[1], rate = obsGamma$estimate[2])
                  s[noRain] <- 0
            } else {
                  warning("There is at least one location without rainfall above the threshold.\n In this (these) location(s) no bias correction has been applied.")
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
#' @param precip Logical indicating if o, p, s is precipitation data.
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm). 
#' @param extrapolation Character indicating the extrapolation method to be applied to correct values in  
#' \code{"s"} that are out of the range of \code{"p"}. Extrapolation is applied only to the \code{"eqm"} method, 
#' thus, this argument is ignored if other bias correction method is selected. Default is \code{"none"} (do not extrapolate).
#' @keywords internal
#' @importFrom stats approxfun ecdf quantile
#' @author S. Herrera and M. Iturbide

eqm <- function(o, p, s, precip, pr.threshold, n.quantiles, extrapolation){
      if (precip == TRUE) {
            threshold <- pr.threshold
            if (any(!is.na(o))) {
                  params <-  norain(o, p, threshold)
                  p <- params$p
                  nP <- params$nP
                  Pth <- params$Pth
            } else {
                  nP = NULL
            }
            if (is.null(nP)) {
                  smap <- rep(NA, length(s))
            } else if (any(!is.na(p)) & any(!is.na(o))) {
                  if (length(which(p > Pth)) > 0) { 
                        noRain <- which(s <= Pth & !is.na(s))
                        rain <- which(s > Pth & !is.na(s))
                        drizzle <- which(s > Pth  & s  <= min(p[which(p > Pth)], na.rm = TRUE) & !is.na(s))
                        eFrc <- tryCatch({ecdf(s[rain])}, error = function(err) {stop("There are not precipitation days in newdata for the step length selected in one or more locations. Try to enlarge the window step")})
                        if (length(rain) > 0) {
                              if (is.null(n.quantiles)) n.quantiles <- length(p)
                              bins <- n.quantiles
                              qo <- quantile(o[which(o > threshold & !is.na(o))], prob = seq(1/bins,1 - 1/bins,1/bins), na.rm = T)
                              qp <- quantile(p[which(p > Pth)], prob = seq(1/bins,1 - 1/bins,1/bins), na.rm = T)
                              p2o <- approxfun(qp, qo, method = "linear")
                              smap <- s
                              smap[rain] <- p2o(s[rain])
                              # Linear extrapolation was discarded due to lack of robustness 
                              if (extrapolation == "constant") {
                                    smap[rain][which(s[rain] > max(qp, na.rm = TRUE))] <- s[rain][which(s[rain] > max(qp, na.rm = TRUE))] + (qo[length(qo)] - qp[length(qo)])
                                    smap[rain][which(s[rain] < min(qp, na.rm = TRUE))] <- s[rain][which(s[rain] < min(qp, na.rm = TRUE))] + (qo[1] - qp[1]) 
                              } else {
                                    smap[rain][which(s[rain] > max(qp, na.rm = TRUE))] <- qo[length(qo)]
                                    smap[rain][which(s[rain] < min(qp, na.rm = TRUE))] <- qo[1]
                              }
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
#' @param precip Logical indicating if o, p, s is precipitation data.
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm). 
#' @param theta numeric indicating  upper threshold (and lower for the left tail of the distributions, if needed) 
#' above which precipitation (temperature) values are fitted to a Generalized Pareto Distribution (GPD). 
#' Values below this threshold are fitted to a gamma (normal) distribution. By default, 'theta' is the 95th 
#' percentile (5th percentile for the left tail). Only for \code{"gpqm"} method.
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
                  params <-  norain(o, p, threshold)
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

#' @title norain
#' @description presprocess to bias correct precipitation data
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm). 
#' @importFrom MASS fitdistr
#' @keywords internal
#' @importFrom stats rgamma
#' @author S. Herrera and M. Iturbide

norain <- function(o, p , threshold){
      nPo <- sum(as.double(o <= threshold & !is.na(o)), na.rm = TRUE)
      nPp <- ceiling(length(p) * nPo / length(o))
      if (nPo >= 0 & nPo < length(o)) {
            ix <- sort(p, decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
            Ps <- sort(p, decreasing = FALSE, na.last = NA)
            Pth <- Ps[nPp + 1]
            if (Pth <= threshold) {
                  Os <- sort(o, decreasing = FALSE, na.last = NA)
                  indP <- which(Ps > threshold & !is.na(Ps))
                  if (length(indP) == 0) {
                        indP <- max(which(!is.na(Ps)))
                        indO <- min(c(length(Os), ceiling(length(Os) * indP/length(Ps))))
                  } else {
                        indP <- min(which(Ps > threshold & !is.na(Ps)))
                        indO <- ceiling(length(Os) * indP/length(Ps))
                  }
                  # [Shape parameter Scale parameter]
                  if (length(unique(Os[(nPo + 1):indO])) < 6) {
                        Ps[(nPp + 1):indP] <- mean(Os[(nPo + 1):indO], na.rm = TRUE)
                  } else {
                        auxOs <- Os[(nPo + 1):indO]
                        auxOs <- auxOs[which(!is.na(auxOs))]
                        auxGamma <- fitdistr(auxOs, "gamma")
                        Ps[(nPp + 1):indP] <- rgamma(indP - nPp, auxGamma$estimate[1], rate = auxGamma$estimate[2])
                  }
                  Ps <- sort(Ps, decreasing = FALSE, na.last = NA)
            }
            if (nPo > 0) {
                  ind <- min(nPp, length(p))
                  Ps[1:ind] <- 0
            }
            p[ix] <- Ps
      } else {
            if (nPo == length(o)) {
                  ix <- sort(p, decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
                  Ps <- sort(p, decreasing = FALSE, na.last = NA)
                  Pth <- Ps[nPp]
                  ind <- min(nPp, length(p))
                  Ps[1:ind] <- 0
                  p[ix] <- Ps
            }
      }
      return(list("nP" = c(nPo,nPp), "Pth" = Pth, "p" = p)) 
}

#end

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
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation.
#' @author B. Szabo-Takacs

loci <- function(o, p, s, precip, pr.threshold){
      if (precip == FALSE) { 
            stop("method loci is only applied to precipitation data")
      } else {
            threshold <- pr.threshold
            l <- length(which(o > threshold))
            gcmr <- rev(sort(s))
            Pgcm <- gcmr[l + 1]
            # local scaling factor
            mobs <- mean(o[which(o > threshold)], na.rm = TRUE)
            mgcm <- mean(s[which(s > Pgcm)], na.rm = TRUE)
            scaling <- (mobs - threshold) / (mgcm - Pgcm)
            GCM <- (scaling*(s - Pgcm)) + threshold
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
#' @importFrom transformeR subsetGrid redim getShape bindGrid.time
#' @author J Bedia

flatMemberDim <- function(grid) {
      grid <- redim(grid, member = TRUE)     
      n.mem.join <- getShape(grid, "member")
      n.time.join <- getShape(grid, "time")
      aux.ltime <- lapply(1:n.mem.join, function(x) {
            subsetGrid(grid, members = x)
      })
      out <- do.call("bindGrid.time", aux.ltime)
      attr(out, "orig.mem.shape") <- n.mem.join
      attr(out, "orig.time.shape") <- n.time.join
      return(out)
}

#' @title Recover member multimember structure
#' @description Recover member multimember structure after application of \code{\link{flatMemberDim}}
#' @param flat.grid A \dQuote{flattened} grid used as predictor in \code{biasCorrection} (the 'pred' object)
#' @param bc.grid The bias-corrected output (the 'bc' object), still without its member structure 
#' @param newdata The 'newdata' object, needed to recover relevant metadata (i.e. initialization dates and member names)
#' @return A (bias-corrected) multimember grid
#' @keywords internal
#' @importFrom transformeR subsetDimension bindGrid.member
#' @seealso \code{\link{flatMemberDim}}, for \dQuote{flattening} the member structure
#' @author J Bedia

recoverMemberDim <- function(flat.grid, bc.grid, newdata) {
      pred <- flat.grid
      bc <- bc.grid
      nmem <- attr(pred, "orig.mem.shape")
      ntimes <- attr(pred, "orig.time.shape")
      bc$Dates <- lapply(bc$Dates, "rep", nmem)
      aux.list <- lapply(1:nmem, function(m) {
            aux <- subsetDimension(grid = bc, dimension = "time", indices = ((m - 1) * ntimes + 1):(m * ntimes))
            aux$InitializationDates <- newdata$InitializationDates[[m]]
            aux$Members <- newdata$Members[[m]]
            return(aux)
      })
      do.call("bindGrid.member", aux.list)
}

