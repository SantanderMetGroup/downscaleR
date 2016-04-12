#' @title Bias correction methods
#' @description Implementation of several standard bias correction methods
#'
#' @template templateObsPredSim
#' @param method method applied. Current accepted values are \code{"eqm"}, \code{"delta"},
#'  \code{"scaling"}, \code{"gqm"} and \code{"gpqm"}. See details.
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm). See details on bias correction for precipitation.
#' @param window vector of length = 2 specifying the time window width used to calibrate and the target days (days that are being corrected).
#'  The window is centered on the target day/s. 
#' Default to \code{NULL}, which considers the whole period available.
#' 
#' @param scaling.type Character indicating the type of the scaling method. Options are \code{"additive"} (default)
#' or \code{"multiplicative"} (see details). This argument is ignored if \code{"scaling"} is not 
#' selected as the bias correction method.
#' 
#' @param extrapolation Character indicating the extrapolation method to be applied to correct values in  
#' \code{"sim"} that are out of the range of \code{"x"}. Extrapolation is applied only to the \code{"eqm"} method, 
#' thus, this argument is ignored if other bias correction method is selected. Default is \code{"no"} (do not extrapolate).
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
#'  \code{"pr"}, \code{"tp"} (this is the standard name defined in the vocabulary (\code{\link[loadeR]{showVocabulary}}), \code{"precipitation"} or \code{"precip"}.
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
#' eqm1 <- biasCorrection(y = y, x = x, newdata = x, method = "eqm", extrapolation = "no", window = NULL)
#' eqm1win <- biasCorrection(y = y, x = x, newdata = x, method = "eqm", extrapolation = "no", window = c(90, 1))
#' NCEP_Igueldo <- subsetGrid(x, latLim = y$xyCoords[,2], lonLim = y$xyCoords[,1])
#' 
#' par(mfrow = c(1,3))
#' qqplot(y$Data, NCEP_Igueldo$Data)
#' lines(c(0,100), c(0,100))
#' qqplot(y$Data, eqm1$Data)
#' lines(c(0,100), c(0,100))
#' qqplot(y$Data, eqm1win$Data)
#' lines(c(0,100), c(0,100))
#' }

biasCorrection <- function(y, x, newdata, method = c("delta", "scaling", "eqm", "gqm", "gpqm"),
                           window = NULL,
                           scaling.type = c("additive", "multiplicative"),
                           pr.threshold = 1, extrapolation = c("no", "constant"), theta = .95){
      obs <- y
      pred <- x
      sim <- newdata
      if (!any(grepl(obs$Variable$varName,c("pr","tp","precipitation","precip")))){
            precip <- FALSE    
      }else{
            precip <- TRUE
      } 
      
      if("station" %in% attr(obs$Data, "dimensions")){
            station <- TRUE
            obs <- redim(obs)
            x <- obs$xyCoords[,1]
            y <- obs$xyCoords[,2]
            ind <- cbind(1:dim(obs$Data)[2], rep(1, dim(obs$Data)[2]), 1:dim(obs$Data)[2])
      }else{
            station <- FALSE
            x <- obs$xyCoords$x
            y <- obs$xyCoords$y
            ind1 <- expand.grid(1:length(y), 1:length(x))
            ind <- cbind(ind1, ind1[,2])
      }
      
      bc <- obs
      pred <- redim(pred)
      sim <- redim(sim)
      
            n.run <- dim(sim$Data)[1]
            n.mem <- dim(sim$Data)[2]
            
            run <- array(dim = c(1, n.mem, dim(sim$Data)[-(1:2)][1], dim(obs$Data)[-1]))
            lrun <- lapply(1:n.run, function(k){    # loop for runtimes
                 
                  pre <- subsetGrid(pred, runtime = k, drop = FALSE)
                  si <- subsetGrid(sim, runtime = k, drop = FALSE)
                  
#                   if(multi.member == TRUE){
                        mem <- array(dim = c(1, 1, dim(sim$Data)[-(1:2)][1], dim(obs$Data)[-1]))
                        lmem <- lapply(1:n.mem, function(l){ # loop for members
                        
                              p <-subsetGrid(pred, members = l, drop = FALSE)
                              s <-subsetGrid(sim, members = l, drop = FALSE)
                              message("[", Sys.time(), "] Bias correcting member ", l, " out of ", n.mem, ".")
                              if(is.null(window)){
                        # Apply bias correction methods
                                          for(i in 1:nrow(ind)){
                              
                                                mem[,,,ind[i,1],ind[i,2]] <- biasCorrection1D(obs$Data[,ind[i,1],ind[i,2]], 
                                                            subsetGrid(p, latLim = y[ind[i,1]], lonLim = x[ind[i,3]], outside = T, drop = FALSE)$Data[1,1,,1,1], 
                                                            subsetGrid(s, latLim = y[ind[i,1]], lonLim = x[ind[i,3]], outside = T, drop = FALSE)$Data[1,1,,1,1],
                                                            method = method,
                                                            scaling.type = scaling.type,
                                                            precip = precip,
                                                            pr.threshold = pr.threshold,
                                                            extrapolation = extrapolation,
                                                            theta = theta)
                                          }                                                                             
                              }else{
                                    step <- window[2]
                                    window <- window[1]+step
                                    message("Correcting windows")
                                    datesList <- as.POSIXct(obs$Dates$start, tz = "GMT", format = "%Y-%m-%d")
                                    yearList <- unlist(strsplit(as.character(datesList), "[-]"))
                                    dayListObs <- array(data = c(as.numeric(yearList[seq(2,length(yearList),3)]),as.numeric(yearList[seq(3,length(yearList),3)])), dim = c(length(datesList),2))
                                    dayList <- unique(dayListObs,index.return = datesList)
                                    indDays <- array(data = NaN, dim = c(length(datesList),1))
                                    for (d in 1:dim(dayList)[1]){
                                          indDays[which(sqrt((dayListObs[,1]-dayList[d,1])^2+(dayListObs[,2]-dayList[d,2])^2)==0)] <- d
                                    }
                                    datesList <- as.POSIXct(sim$Dates$start, tz="GMT", format="%Y-%m-%d")
                                    yearList<-unlist(strsplit(as.character(datesList), "[-]"))
                                    dayListSim <- array(data = c(as.numeric(yearList[seq(2,length(yearList),3)]),as.numeric(yearList[seq(3,length(yearList),3)])), dim = c(length(datesList),2))
                                    indDaysSim <- array(data = NaN, dim = c(length(datesList),1))
                                    for (d in 1:dim(dayList)[1]){
                                          indDaysSim[which(sqrt((dayListSim[,1]-dayList[d,1])^2+(dayListSim[,2]-dayList[d,2])^2)==0)] <- d
                                    }
                                    for(i in 1:nrow(ind)){
                                          end <- indDaysSim
                                          steps <- floor(dim(dayList)[1]/step)
                                          for (j in 1:steps){
                                                days <- ((j-1)*step+1):((j-1)*step+step)
                                                if(j == steps) days <- days[1]:dim(dayList)[1]
                                                indObs <- lapply(1:length(days), function(h){
                                                      which(indDays == days[h])
                                                })
                                                indObs <- sort(do.call("abind", indObs))
                                                indObsWindow <- array(data = NA, dim = c((window)*length(indObs)/step,1))
                                                breaks <- c(which(diff(indObs) != 1), length(indObs))
                                                for (d in 1:length(breaks)){
                                                            if(d == 1) {
                                                                  piece <- indObs[1:breaks[1]]
                                                            }else{
                                                                  piece <- indObs[(breaks[d-1]+1): breaks[d]]
                                                            }
                                                            suppressWarnings(indObsWindow[((d-1)*window+1):(d*window)] <- (min(piece, na.rm = T)-floor((window-step)/2)):(max(piece, na.rm = T)+floor((window-step)/2)))
                                                }
                                                indObsWindow[which(indObsWindow <= 0)] <- 1
                                                indObsWindow[which(indObsWindow >  length(indDays))] <- length(indDays)
                                                indObsWindow <- unique(indObsWindow)
                                                indSim <- lapply(1:length(days), function(h){
                                                      which(indDaysSim == days[h])
                                                })
                                                indSim <- sort(do.call("abind", indSim))
                                          # Apply bias correction methods
                                          
                                                o1 <- obs$Data[indObsWindow,ind[i,1],ind[i,2]]
                                                p1 <- subsetGrid(p, latLim = y[ind[i,1]], lonLim = x[ind[i,3]], outside = T, drop = FALSE)$Data[1,1,,1,1][indObsWindow]
                                                s1 <- subsetGrid(s, latLim = y[ind[i,1]], lonLim = x[ind[i,3]], outside = T, drop = FALSE)$Data[1,1,,1,1][indSim]
                                                end[indSim,] <- biasCorrection1D(o1, p1, s1, 
                                                                 method = method,
                                                                 scaling.type = scaling.type,
                                                                 precip = precip,
                                                                 pr.threshold = pr.threshold,
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
            ##########################
            bc <- redim(bc, drop = TRUE)
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
#' thus, this argument is ignored if other bias correction method is selected. Default is \code{"no"} (do not extrapolate).
#' @param theta numeric indicating  upper threshold (and lower for the left tail of the distributions, if needed) 
#' above which precipitation (temperature) values are fitted to a Generalized Pareto Distribution (GPD). 
#' Values below this threshold are fitted to a gamma (normal) distribution. By default, 'theta' is the 95th 
#' percentile (5th percentile for the left tail). Only for \code{"gpqm"} method.
#' @keywords internal
#' @author M. Iturbide

biasCorrection1D <- function(o, p, s, method = c("delta", "scaling", "eqm", "gqm", "gpqm"), 
                             scaling.type = c("additive", "multiplicative"), precip = NULL, 
                             pr.threshold = 1, extrapolation = c("no", "constant"), theta = .95){
      if(method == "delta"){
            bc1d <- delta(o, p, s)
      }else if(method == "scaling"){
            bc1d <- scaling(o, p, s, scaling.type)
      }else if(method == "eqm"){
            bc1d <- eqm(o, p, s, precip, pr.threshold, extrapolation)
      }else if(method == "gqm"){
            bc1d <- gqm(o, p, s, precip, pr.threshold)
      }else if(method == "gpqm"){
            bc1d <- gpqm(o, p, s, precip, pr.threshold, theta)
      }
      
     return(bc1d) 
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
      if (scaling.type == "additive"){
          corrected <- s - mean(p) + mean(o)
            
      }else if(scaling.type == "multiplicative"){
           corrected <- (s/mean(p))* mean(o)
            
      }
      return(corrected)
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
#' @keywords internal
#' @author S. Herrera and M. Iturbide
gqm <- function(o, p, s, precip, pr.threshold){
                
      if (precip == FALSE) {
            stop("method gqm is only applied to precipitation data")
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
           if (is.null(nP)){
                 s <- o
           }else if(nP < length(o)){
                  ind <- which(o > threshold & !is.na(o))
                  obsGamma <- fitdistr(o[ind],"gamma")
                  ind <- which(p > 0 & !is.na(p))
                  prdGamma <- fitdistr(p[ind],"gamma")
                  rain <- which(s > Pth & !is.na(s))
                  noRain <- which(s <= Pth & !is.na(s))
                  auxF <- pgamma(s[rain], prdGamma$estimate[1], rate = prdGamma$estimate[2])
                  s[rain] <- qgamma(auxF, obsGamma$estimate[1], rate = obsGamma$estimate[2])
                  s[noRain] <- 0
                  warningNoRain <- FALSE
           }else{
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
#' thus, this argument is ignored if other bias correction method is selected. Default is \code{"no"} (do not extrapolate).
#' @keywords internal
#' @author S. Herrera and M. Iturbide
eqm <- function(o, p, s, precip, pr.threshold, extrapolation){
      if (precip == TRUE) {
            threshold <- pr.threshold
                        if (any(!is.na(o))){
                              params <-  norain(o, p, threshold)
                              p <- params$p
                              nP <- params$nP
                              Pth <- params$Pth
                        }else{
                              nP = NULL
                        }
           ################################################################################################-
            if(is.null(nP)){
                  s <- o
            }else if(any(!is.na(p)) & any(!is.na(o))){
                  if (length(which(p > Pth)) > 0){ 
                        ##
                        ePrd<-ecdf(p[which(p > Pth)]) 
                        ##
                        noRain <- which(s<=Pth & !is.na(s))
                        rain <- which(s > Pth & !is.na(s))
                        drizzle <- which(s > Pth  & s  <= min(p[which(p > Pth)], na.rm = TRUE) & !is.na(s))
                        if (length(rain) > 0){
                              eFrc <- ecdf(s[rain]) #eFrc??
                              ##############################################################################-
                              if(extrapolation == "constant"){
                                    exupsim <- which(s[rain] > max(range(p, na.rm = TRUE))) #
                                    exdwnsim <- which(s[rain] < min(range(p, na.rm = TRUE)))  
                                    ex <- c(exupsim,exdwnsim)
                                    exC <- setdiff(1:length(rain),ex)
                                    if (length(exupsim)>0){
                                          extmin <- max(p, na.rm = TRUE)
                                          dif <- extmin - max(o[which(o > threshold & !is.na(o))], na.rm = TRUE)
                                          s[rain[exupsim]] <- s[rain[exupsim]] - dif
                                    }
                                    if (length(exdwnsim) > 0){
                                          extmin <- min(p, na.rm = TRUE)
                                          dif <- extmin - min(o[which(o > threshold & !is.na(obs))], na.rm = TRUE)
                                          s[rain[exdwnsim]] <- s[rain[exdwnsim]] - dif
                                    }
                                    s[rain[exC]] <- quantile(o[which(o > threshold & !is.na(o))], probs = ePrd(s[rain[exC]]), na.rm = TRUE, type = 4)
                              }else{
                                    s[rain]<-quantile(o[which(o > threshold & !is.na(o))], probs = ePrd(s[rain]), na.rm = TRUE, type = 4)
                              }
                              ####################################################################################-
                        }
                        if (length(drizzle)>0){
                              s[drizzle]<-quantile(s[which(s > min(p[which(p > Pth)], na.rm = TRUE) & !is.na(s))], probs = eFrc(s[drizzle]), na.rm = TRUE, type = 4)
                        }
                        s[noRain]<-0
                  }else{
                        noRain<-which(s <= Pth & !is.na(s))
                        rain<-which(s > Pth & !is.na(s))
                        if (length(rain)>0){
                              eFrc<-ecdf(s[rain])
                              s[rain]<-quantile(o[which(o > threshold & !is.na(o))], probs = eFrc(s[rain]), na.rm = TRUE, type = 4)
                        }
                        s[noRain]<-0
                  }
            }
       }else{
             if(all(is.na(o))){
                   s <- o
             }else if (any(!is.na(p)) & any(!is.na(o))) {
                        ePrd <- ecdf(p)
                              #################-
                              if(extrapolation == "constant"){
                                    exupsim <- which(s > max(range(p, na.rm = TRUE))) 
                                    exdwnsim <- which(s < min(range(p, na.rm = TRUE)))  
                                    ex <- c(exupsim,exdwnsim)
                                    exC <- setdiff(1:length(s),ex)
                                    if (length(exupsim)>0){
                                          extmin <- max(p, na.rm = TRUE)
                                          dif <- extmin - max(o, na.rm = TRUE)
                                          s[exupsim] <- s[exupsim] - dif
                                    }
                                    if (length(exdwnsim)>0){
                                          extmin <- min(p, na.rm = TRUE)
                                          dif <- extmin - min(o, na.rm = TRUE)
                                          s[exdwnsim] <- s[exdwnsim] - dif
                                    }
                                    s[exC] <- quantile(o, probs = ePrd(s[exC]), na.rm = TRUE, type = 4)
                              }else{
                                    s <- quantile(o, probs = ePrd(s), na.rm = TRUE, type = 4)
                              }
                        
                              #################-
            }
      }
      return(s)
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
#' @importFrom evd qgpd
#' @importFrom evd pgpd
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
                     s <- o
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
