#' @title Bias correction methods
#' @description Implementation of several standard bias correction methods
#'
#' @template templateObsPredSim
#' @param method method applied. Current accepted values are \code{"eqm"}, \code{"delta"},
#'  \code{"scaling"}, \code{"gqm"} and \code{"gpqm"}. See details.
#' @param precipitation Logical indicating if the data to be corrected is precipitation data.
#' @param cross.val Should cross-validation be performed? methods available are leave-one-out ("loocv") and k-fold ("kfold"). The default 
#' option ("none") does not perform cross-validation.
#' @param folds Only requiered if cross.val = "kfold". A list of vectors, each containing the years to be grouped in 
#' the corresponding fold.
#' @param wet.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm). See details on bias correction for precipitation.
#' @param window vector of length = 2 specifying the time window width used to calibrate and the target days (days that are being corrected).
#'  The window is centered on the target day/s. 
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
#' 
#' @param theta numeric indicating  upper threshold (and lower for the left tail of the distributions, if needed) 
#' above which precipitation (temperature) values are fitted to a Generalized Pareto Distribution (GPD). 
#' Values below this threshold are fitted to a gamma (normal) distribution. By default, 'theta' is the 95th 
#' percentile (5th percentile for the left tail). Only for \code{"gpqm"} method.
#'  
#' @details
#' 
#' The methods available are \code{"eqm"}, \code{"delta"}, 
#' \code{"scaling"}, \code{"gqm"}, \code{"gpqm"} (the two latter used only for precipitation).
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
#' @seealso \code{\link{isimip}} for a trend-preserving method of model calibration
#' @return A calibrated grid of the same spatio-temporal extent than the input \code{"y"}
#' @family downscaling
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
#' @author S. Herrera and M. Iturbide
#' @export
#' @examples \dontrun{
#' data(VALUE_Igueldo_tp)
#' data(NCEP_Iberia_tp)
#' y <- VALUE_Igueldo_tp
#' x <- NCEP_Iberia_tp
#' 
#' eqm1 <- biasCorrection(y = y, x = x, newdata = x,
#'                        method = "eqm",
#'                        extrapolation = "none",
#'                        window = NULL)
#' eqm1win <- biasCorrection(y = y, x = x, newdata = x,
#'                           method = "eqm",
#'                           extrapolation = "none",
#'                           window = c(90, 10))
#' NCEP_Igueldo <- subsetGrid(x, latLim = y$xyCoords[,2], lonLim = y$xyCoords[,1])
#' par(mfrow = c(1,3))
#' qqplot(y$Data, NCEP_Igueldo$Data)
#' lines(c(0,100), c(0,100))
#' qqplot(y$Data, eqm1$Data)
#' lines(c(0,100), c(0,100))
#' qqplot(y$Data, eqm1win$Data)
#' lines(c(0,100), c(0,100))
#' par(mfrow = c(1,1))
#' }


biasCorrection <- function(y, x, newdata, precipitation = FALSE,
                           method = c("delta", "scaling", "eqm", "gqm", "gpqm"),
                           cross.val = c("none", "loocv", "kfold"),
                           folds = NULL,
                           window = NULL,
                           scaling.type = c("additive", "multiplicative"),
                           wet.threshold = 1, n.quantiles = NULL, extrapolation = c("none", "constant"), 
                           theta = .95){
      method <- match.arg(method, choices = c("delta", "scaling", "eqm", "gqm", "gpqm"))
      cross.val <- match.arg(cross.val, choices = c("none", "loocv", "kfold"))
      scaling.type <- match.arg(scaling.type, choices = c("additive", "multiplicative"))
      extrapolation <- match.arg(extrapolation, choices = c("none", "constant"))
      if(cross.val == "none"){
            output <- biasCorrectionXD(y = y, x = x, newdata = newdata, 
                                       precipitation = precipitation,
                                       method = method,
                                       window = window,
                                       scaling.type = scaling.type,
                                       pr.threshold = wet.threshold, 
                                       n.quantiles = n.quantiles, 
                                       extrapolation = extrapolation, 
                                       theta = theta)
      }else{
            if (!is.null(newdata)) {
                  message("'newdata' will be ignored for cross-validation")
            }
            if(cross.val == "loocv"){
                  years <- as.list(unique(getYearsAsINDEX(x)))
            }else if(cross.val == "kfold" & !is.null(folds)){
                  years <- folds
            }else if(cross.val == "kfold" & is.null(folds)){
                  stop("Please, specify folds for kfold cross validation")
            }
            output.list <- lapply(1:length(years), function(i) {
                        target.year <-years[[i]]
                        rest.years <- setdiff(unlist(years), target.year)
                        yy <- redim(y, member = F)
                        yy <- subsetGrid(yy, years = rest.years, drop = FALSE)
                        yy <- redim(yy, drop = T)
                        if(length(getDim(yy)) == 1){attr(yy$Data, "dimensions") <- c(getDim(yy), "station")}
                        newdata2 <- subsetGrid(x, years = target.year)
                        xx <- subsetGrid(x, years = rest.years)
                        message("Validation ", i, ", ", length(unique(years)) - i, " remaining")
                        biasCorrectionXD(y = yy, x = xx, newdata = newdata2, precipitation = precipitation,
                                          method = method,
                                          window = window,
                                          scaling.type = scaling.type,
                                          pr.threshold = wet.threshold, n.quantiles = n.quantiles, extrapolation = extrapolation, 
                                          theta = theta)
                  })
                  al <- which(getDim(x) == "time")
                  Data <- sapply(output.list, function(n) unname(n$Data), simplify = F)
                  bindata <- unname(do.call("abind", c(Data, along = al)))
                  output <- output.list[[1]]
                  output$Data <- bindata
                  output$Dates <- x$Dates
      }
      return(output)
}

biasCorrectionXD <- function(y, x, newdata, precipitation, 
                           method = c("delta", "scaling", "eqm", "gqm", "gpqm"),
                           window = NULL,
                           scaling.type = c("additive", "multiplicative"),
                           pr.threshold = 1, n.quantiles = NULL, extrapolation = c("none", "constant"), 
                           theta = .95){
      method <- match.arg(method, choices = c("delta", "scaling", "eqm", "gqm", "gpqm"))
      scaling.type <- match.arg(scaling.type, choices = c("additive", "multiplicative"))
      extrapolation <- match.arg(extrapolation, choices = c("none", "constant"))
      obs <- y
      pred <- x
      sim <- newdata
#       if (!any(grepl(obs$Variable$varName,c("pr","tp","precipitation","precip")))) {
#             precip <- FALSE    
#         } else {
#             precip <- TRUE
#         }
      precip <- precipitation
      if ("station" %in% attr(obs$Data, "dimensions")) {
            station <- TRUE
            obs <- redim(obs, member = F)
            x <- obs$xyCoords[,1]
            y <- obs$xyCoords[,2]
            ito <- which(getDim(obs) == "time")
            ind <- cbind(1:dim(obs$Data)[2], rep(1, dim(obs$Data)[2]), 1:dim(obs$Data)[2])
      } else {
            station <- FALSE
            x <- obs$xyCoords$x
            y <- obs$xyCoords$y
            ind1 <- expand.grid(1:length(y), 1:length(x))
            ind <- cbind(ind1, ind1[,2])
      }
      
      bc <- obs
      pred <- redim(pred, member = T, runtime = T)
      sim <- redim(sim, member = T, runtime = T)
      itp <- which(getDim(pred) == "time")
      ito <- which(getDim(obs) == "time")
      if (dim(obs$Data)[ito] != dim(pred$Data)[itp]) stop("y and x do not have the same time series length")
      
      n.run <- dim(sim$Data)[1]
      n.mem <- dim(sim$Data)[2]
      
      if (method == "delta") {
            run <- array(dim = c(1, n.mem, dim(obs$Data)))
      } else {
            its <- which(getDim(sim) == "time")
            run <- array(dim = c(1, n.mem, dim(sim$Data)[its], dim(obs$Data)[-ito]))
      }
      lrun <- lapply(1:n.run, function(k) {    # loop for runtimes
            
            pre <- subsetGrid(pred, runtime = k, drop = FALSE)
            si <- subsetGrid(sim, runtime = k, drop = FALSE)
            
            #                   if(multi.member == TRUE){
            if (method == "delta") {
                  mem <- array(dim = c(1, 1, dim(obs$Data)))
            } else {
                  mem <- array(dim = c(1, 1, dim(sim$Data)[its], dim(obs$Data)[-ito]))
            }
            lmem <- lapply(1:n.mem, function(l){ # loop for members
                  
                  p <- subsetGrid(pred, members = l, drop = FALSE)
                  s <- subsetGrid(sim, members = l, drop = FALSE)
                  message("[", Sys.time(), "] Bias correcting member ", l, " out of ", n.mem, ".")
                  if (is.null(window)) {
                        # Apply bias correction methods
                        for (i in 1:nrow(ind)) {
                              suppressWarnings(
                                   mem[,,,ind[i,1],ind[i,2]] <- biasCorrection1D(obs$Data[,ind[i,1],ind[i,2]],
                                                                            subsetGrid(p, latLim = y[ind[i,1]], lonLim = x[ind[i,3]], outside = T, drop = FALSE)$Data[1,1,,1,1],
                                                                            subsetGrid(s, latLim = y[ind[i,1]], lonLim = x[ind[i,3]], outside = T, drop = FALSE)$Data[1,1,,1,1], 
                                                                            method = method,
                                                                            scaling.type = scaling.type,
                                                                            precip = precip,
                                                                            pr.threshold = pr.threshold,
                                                                            n.quantiles = n.quantiles,
                                                                            extrapolation = extrapolation,
                                                                            theta = theta)
                              )
                        }                                                                             
                  } else {
                        step <- window[2]
                        window <- window[1] + step
                        message("Correcting windows")
                        datesList <- as.POSIXct(obs$Dates$start, tz = "GMT", format = "%Y-%m-%d")
                        yearList <- unlist(strsplit(as.character(datesList), "[-]"))
                        dayListObs <- array(data = c(as.numeric(yearList[seq(2,length(yearList),3)]),as.numeric(yearList[seq(3,length(yearList),3)])), dim = c(length(datesList),2))
                        dayList <- unique(dayListObs,index.return = datesList)
                        indDays <- array(data = NaN, dim = c(length(datesList),1))
                        for (d in 1:dim(dayList)[1]) {
                              indDays[which(sqrt((dayListObs[,1] - dayList[d,1]) ^ 2 + (dayListObs[,2] - dayList[d,2]) ^ 2) == 0)] <- d
                        }
                        datesList <- as.POSIXct(sim$Dates$start, tz = "GMT", format = "%Y-%m-%d")
                        yearList <- unlist(strsplit(as.character(datesList), "[-]"))
                        dayListSim <- array(data = c(as.numeric(yearList[seq(2,length(yearList),3)]),as.numeric(yearList[seq(3,length(yearList),3)])), dim = c(length(datesList),2))
                        indDaysSim <- array(data = NaN, dim = c(length(datesList),1))
                        for (d in 1:dim(dayList)[1]) {
                              indDaysSim[which(sqrt((dayListSim[,1] - dayList[d,1]) ^ 2 + (dayListSim[,2] - dayList[d,2]) ^ 2) == 0)] <- d
                        }
                        for (i in 1:nrow(ind)) {
                              if (method == "delta") {
                                    end <- indDays
                              } else {
                                    end <- indDaysSim
                              }
                              steps <- floor(dim(dayList)[1]/step)
                              for (j in 1:steps) {
                                    days <- ((j - 1) * step + 1):((j - 1) * step + step)
                                    if (j == steps) days <- days[1]:dim(dayList)[1]
                                    indObs <- lapply(1:length(days), function(h){
                                          which(indDays == days[h])
                                    })
                                    indObs <- sort(do.call("abind", indObs))
                                    indObsWindow <- array(data = NA, dim = c((window)*length(indObs)/step,1))
                                    breaks <- c(which(diff(indObs) != 1), length(indObs))
                                    for (d in 1:length(breaks)) {
                                          if (d == 1) {
                                                piece <- indObs[1:breaks[1]]
                                          } else {
                                                piece <- indObs[(breaks[d - 1] + 1):breaks[d]]
                                          }
                                          suppressWarnings(indObsWindow[((d - 1)*window + 1):(d*window)] <- (min(piece, na.rm = T) - floor((window - step) / 2)):(max(piece, na.rm = T) + floor((window - step) / 2)))
                                    }
                                    indObsWindow[which(indObsWindow <= 0)] <- 1
                                    indObsWindow[which(indObsWindow >  length(indDays))] <- length(indDays)
                                    indObsWindow <- unique(indObsWindow)
                                    indSim <- lapply(1:length(days), function(h){
                                          which(indDaysSim == days[h])
                                    })
                                    indSim <- sort(do.call("abind", indSim))
                                    # Apply bias correction methods
                                    
                                    suppressWarnings(
                                          if (method == "delta") {
                                                o1 <- obs$Data[indObs,ind[i,1],ind[i,2]]
                                          } else {
                                                o1 <- obs$Data[indObsWindow,ind[i,1],ind[i,2]]
                                          }
                                    )
                                    suppressWarnings(
                                          p1 <- subsetGrid(p, latLim = y[ind[i,1]], lonLim = x[ind[i,3]], outside = T, drop = FALSE)$Data[1,1,,1,1][indObsWindow]
                                    )
                                    suppressWarnings(
                                          s1 <- subsetGrid(s, latLim = y[ind[i,1]], lonLim = x[ind[i,3]], outside = T, drop = FALSE)$Data[1,1,,1,1][indSim]
                                    )
                                    if (method == "delta") indSim <- indObs
                                    end[indSim,] <- biasCorrection1D(o1, p1, s1, 
                                                                     method = method,
                                                                     scaling.type = scaling.type,
                                                                     precip = precip,
                                                                     pr.threshold = pr.threshold,
                                                                     n.quantiles = n.quantiles,
                                                                     extrapolation = extrapolation,
                                                                     theta = theta)
                              }
                              mem[,,,ind[i,1],ind[i,2]] <- end
                        }
                  }
                  return(mem)
            }) #end loop for members
            #                   }else{
            #                         ############# members jointly
            #                   }
            run[,,,,] <- abind(lmem, along = 2)
            return(run)
      }) #end loop for runtimes
      bc$Data <- unname(abind(lrun, along = 1))
      attr(bc$Data, "dimensions") <- attr(sim$Data, "dimensions")
      attr(bc$Variable, "correction") <- method
      ##########################
      bc <- redim(bc, runtime = TRUE, drop = TRUE)
      message("[", Sys.time(), "] Done.")
      ##########################
      bc$Dates <- sim$Dates
      bc$InitializationDates <- sim$InitializationDates
      bc$Members <- sim$Members
      return(bc)
}
#end


#' @title Bias correction methods on 1D data
#' @description Implementation of several standard bias correction methods
#' @param o A vector (e.g. station data) containing the observed climate data for the training period
#' @param p A vector containing the simulated climate by the model for the training period. 
#' @param s A vector containing the simulated climate for the variable used in \code{p}, but considering the test period.
#' @param method method applied. Current accepted values are \code{"eqm"}, \code{"delta"},
#'  \code{"scaling"}, \code{"gqm"} and \code{"gpqm"}. 
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
      }
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
#end

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
#end


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
            } else if (nP < length(o)) {
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
                  warning("There is at least one location without rainfall above the threshold.\n In this (these) location(s) none bias correction has been applied.")
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
                  } else {## For dry series
                        smap <- s
                        warning('No rainy days in the prediction. Bias correction is not applied') 
#                         noRain<-which(s <= Pth & !is.na(s))
#                         rain<-which(s > Pth & !is.na(s))
#                         smap <- s
#                         if (length(rain)>0){
#                               eFrc<-ecdf(s[rain])
#                               smap[rain]<-quantile(o[which(o > threshold & !is.na(o))], probs = eFrc(s[rain]), na.rm = TRUE, type = 4)
# 
#                         }
#                         smap[noRain]<-0
                  }
            }
      } else {
            if (all(is.na(o))) {
                  smap <- rep(NA, length(s))
            }else if (all(is.na(p))){
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

gpqm <- function(o, p, s, precip, pr.threshold, theta){ 
      if(precip == FALSE){
            stop("method gpqm is only applied to precipitation data")
      }else{
            
            threshold <- pr.threshold
            if (any(!is.na(o))){
                  params <-  norain(o, p, threshold)
                  p <- params$p
                  nP <- params$nP
                  Pth <- params$Pth
            }else{
                  nP = NULL
            }
            if(is.null(nP)){
                  s <- rep(NA, length(s))
            }else if(nP < length(o)){
                  ind <- which(o > threshold & !is.na(o))
                  indgamma <- ind[which(o[ind] < quantile(o[ind], theta))]
                  indpareto <- ind[which(o[ind] >= quantile(o[ind], theta))]
                  obsGQM <- fitdistr(o[indgamma],"gamma")
                  obsGQM2 <- fpot(o[indpareto], quantile(o[ind], theta), "gpd", std.err = FALSE)
                  ind <- which(p > 0 & !is.na(p))
                  indgamma <- ind[which(p[ind] < quantile(p[ind],theta))]
                  indpareto <-ind[which(p[ind] >= quantile(p[ind], theta))]
                  prdGQM <- fitdistr(p[indgamma], "gamma")
                  prdGQM2 <- fpot(p[indpareto], quantile(p[ind], theta), "gpd", std.err = FALSE)
                  rain <- which(s > Pth & !is.na(s))
                  noRain <- which(s <= Pth & !is.na(s))
                  indgammasim <- rain[which(s[rain] < quantile(p[ind], theta))]
                  indparetosim <- rain[which(s[rain] >= quantile(p[ind], theta))]
                  auxF <- pgamma(s[indgammasim], prdGQM$estimate[1], rate = prdGQM$estimate[2])
                  auxF2 <- pgpd(s[indparetosim], loc = 0, scale = prdGQM2$estimate[1], shape = prdGQM2$estimate[2])
                  s[indgammasim] <- qgamma(auxF, obsGQM$estimate[1], rate = obsGQM$estimate[2])
                  s[indparetosim[which(auxF2<1)]] <- qgpd(auxF2[which(auxF2 < 1)], loc = 0, scale = obsGQM2$estimate[1], shape = obsGQM2$estimate[2])
                  s[indparetosim[which(auxF2==1)]] <- max(o[indpareto], na.rm = TRUE)
                  s[noRain] <- 0
                  warningNoRain <- FALSE
            }else{
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
      nP <- sum(as.double(o <= threshold & !is.na(o)), na.rm = TRUE)
      if (nP >= 0 & nP < length(o)){
            ix <- sort(p, decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
            Ps <- sort(p, decreasing = FALSE, na.last = NA)
            Pth <- Ps[nP + 1]
            if (Ps[nP + 1]<=threshold){
                  Os <- sort(o, decreasing = FALSE, na.last = NA)
                  ind <- which(Ps > threshold & !is.na(Ps))
                  if (length(ind)==0){
                        ind <- max(which(!is.na(Ps)))
                        ind <- min(c(length(Os),ind))
                  }else{
                        ind <- min(which(Ps > threshold & !is.na(Ps)))
                  }
                  # [Shape parameter Scale parameter]
                  if (length(unique(Os[(nP + 1):ind])) < 6){
                        Ps[(nP + 1):ind] <- mean(Os[(nP + 1):ind], na.rm = TRUE)
                  }else{
                        auxOs <- Os[(nP + 1):ind]
                        auxOs <- auxOs[which(!is.na(auxOs))]
                        auxGamma <- fitdistr(auxOs,"gamma")
                        Ps[(nP + 1):ind]<-rgamma(ind - nP, auxGamma$estimate[1], rate = auxGamma$estimate[2])
                  }
                  Ps<-sort(Ps, decreasing = FALSE, na.last = NA)
            }
            if (nP > 0){
                  ind <- min(nP, length(p))
                  Ps[1:ind]<-0
            }
            p[ix] <- Ps
            
      }else{
            if (nP == length(o)){
                  ix <- sort(p, decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
                  Ps <- sort(p, decreasing = FALSE, na.last = NA)
                  Pth <- Ps[nP]
                  ind <- min(nP, length(p))
                  Ps[1:ind] <- 0
                  p[ix] <- Ps
            }
      }
      return(list("nP" = nP, "Pth" = Pth, "p" = p)) 
}

#end
