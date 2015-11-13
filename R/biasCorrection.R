
#' @title Bias correction methods
#' @description Implementation of several standard bias correction methods
#' 
#' @template templateObsPredSim
#' @param method method applied. Current accepted values are \code{"eqm"}, \code{"delta"},
#'  \code{"scaling"}, \code{"gqm"} and \code{"gpqm"}. See details.
#' @param multi.member Should members be adjusted sepparately (TRUE, default), or jointly (FALSE)?. 
#' Ignored if the dataset has no members.
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm). See details on bias correction for precipitation.
#' @param window Numeric value specifying the time window width used to calibrate. The window is centered on the target day. 
#' Default to \code{NULL}, which considers the whole period available.
#' 
#' @param scaling.type Character indicating the type of the scaling method. Options are \code{"additive"} (default)
#' or \code{"multiplicative"} (see details). This argument is ignored if \code{"scaling"} is not 
#' selected as the bias correction method.
#' 
#' @param extrapolation Character indicating the extrapolation method to be applied to correct values in  
#' \code{"sim"} that are out of the range of \code{"pred"}. Extrapolation is applied only to the \code{"eqm"} method, 
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
#'\strong{gpqm}
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
#'  \code{"pr"}, \code{"tp"} (this is the standard name defined in the \code{\link{vocabulary}}), \code{"precipitation"} or \code{"precip"}.
#'  Thus, caution must be taken to ensure that the correct bias correction is being undertaken when dealing with
#'  non-standard variables.
#'     
#' 
#' @seealso \code{\link{isimip}} for a trend-preserving method of model calibration
#'  
#' @return A calibrated object of the same spatio-temporal extent of the input field
#'  
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
#' 
#' 
#' @author S. Herrera \email{sixto@@predictia.es}
#' @export
#' @examples \dontrun{
#' # Download VALUE (station data) and NCEP (model data) datasets 
#' dir.create("mydirectory")
#' download.file("http://meteo.unican.es/work/downscaler/data/VALUE_ECA_86_v2.tar.gz", 
#'               destfile = "mydirectory/VALUE_ECA_86_v2.tar.gz")
#' download.file("http://meteo.unican.es/work/downscaler/data/Iberia_NCEP.tar.gz", 
#'               destfile = "mydirectory/Iberia_NCEP.tar.gz")
#' # Extract files from the tar.gz file
#' untar("mydirectory/VALUE_ECA_86_v2.tar.gz", exdir = "mydirectory")
#' untar("mydirectory/NCEP_Iberia.tar.gz", exdir = "mydirectory")
#' # Path to the VALUE dataset and the NCEP ncml file.
#' value <- "mydirectory/VALUE_ECA_86_v2"
#' ncep <- "mydirectory/Iberia_NCEP/Iberia_NCEP.ncml"
#' # Data inventories provides a quick overview of the available data
#' value.inv <- dataInventory(value)
#' ncep.inv <- dataInventory(ncep)
#' str(value.inv)
#' str(ncep.inv)
#' # Load precipitation for boreal winter (DJF) in the train (1991-2000) and test (2001-2010) periods,
#' # for the observations (VALUE) and the Iberia_NCEP datasets
#' obs <- loadStationData(dataset = value,
#'                        var="precip",
#'                        lonLim = c(-12,10), latLim = c(33,47),
#'                        season = c(12,1,2),
#'                        years = 1991:2000)
#' prd <- loadGridData(dataset = ncep,
#'                     var = "tp",
#'                     lonLim = c(-12,10), latLim = c(33,47),
#'                     season = c(12,1,2),
#'                     years = 1991:2000)
#' sim <- loadGridData(dataset = ncep,
#'                     var = "tp",
#'                     lonLim = c(-12,10), latLim = c(33,47),
#'                     season = c(12,1,2),
#'                     years = 2001:2010)
#' # Interpolation of the observations onto the grid of model: we use the method "nearest" and the
#' # 'getGrid' function to ensure spatial consistency:
#' obs <- interpData(obs, new.Coordinates = getGrid(prd), method = "nearest")
#' # Apply the bias correction method:
#' simBC <- biasCorrection (obs, prd, sim, method = "eqm", pr.threshold = 1) # qq-mapping
#' par(mfrow = c(1,2))
#' plotMeanField(sim)
#' plotMeanField(simBC)
#' }

biasCorrection <- function (obs, pred, sim, method = c("eqm", "delta", "scaling", "gqm", "gpqm"), pr.threshold = 1, 
                            multi.member = TRUE, window = NULL, scaling.type = c("additive", "multiplicative"), extrapolation = c("no", "constant"), theta = .95) {
  method <- match.arg(method, choices = c("eqm", "delta", "scaling", "gqm", "gpqm"))
  extrapolation <- match.arg(extrapolation, choices = c("no", "constant"))
  threshold <- pr.threshold
  obj <- getIntersect(obs,pred)
  obs <- obj$obs
  pred <- obj$prd
  if  (length(obj$Dates$start)==0){
    stop("The obs and pred periods do not overlap\nCheck the periods")
  }else{
    if ((length(obj$Dates$start)!=length(obs$Dates$start)) | (length(obj$Dates$start)!=length(pred$Dates$start))){
      warning("The objects obs and pred, used for the method calibration, have different periods.\n Only the overlap period is considered for the calibration.")
    }
  }
  dimObs <- dim(obs$Data)
  obs.time.index <- grep("^time$", attr(obs$Data, "dimensions"))
  dimPred <- dim(pred$Data)
  dimFor <- dim(sim$Data)
  pred.time.index <- grep("^time$", attr(pred$Data, "dimensions"))
  dimDiff <- NULL
  dimPerm <- 1:length(attr(pred$Data, "dimensions"))
  dimPerm[which(dimPerm<=length(dim(obs$Data)))] <- match(attr(obs$Data, "dimensions"), attr(pred$Data, "dimensions"))
  if (length(setdiff(attr(pred$Data, "dimensions"),attr(obs$Data, "dimensions")))>0){
    dimDiff<-setdiff(attr(pred$Data, "dimensions"),attr(obs$Data, "dimensions"))
    for (k in 1:length(dimDiff)) {
      indDiff <- which(grepl(dimDiff[k],attr(pred$Data, "dimensions")))
      dimPerm[length(attr(obs$Data, "dimensions"))+k] <- indDiff
    }
  }
  datesList <- as.POSIXct(obs$Dates$start, tz="GMT", format="%Y-%m-%d")
  yearList<-unlist(strsplit(as.character(datesList), "[-]"))
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
  attrSim <- attr(sim$Data, "dimensions")
  #apply function of calibration
  if (length(dimDiff)==0){
    F <- calibrateProj(obs$Data, aperm(pred$Data,dimPerm), aperm(sim$Data,dimPerm), method = method, varcode = obs$Variable$varName, pr.threshold = threshold, scaling.type = scaling.type, extrapolate = extrapolation, theta = theta)
    if (!any(dimPerm != 1:length(attr(pred$Data, "dimensions")))){
      sim$Data<-F
    }else{
      sim$Data<-aperm(F,c(dimPerm[2:length(dim(obs$Data))],dimPerm[1]))
    }
  }else{
    if ((("member" %in% dimDiff) & isTRUE(multi.member)) | !("member" %in% dimDiff)) {
      indDimDiff <- rep(list(bquote()), length(dimDiff))
      for (d in 1:length(dimDiff)){
        indDimDiff[[d]] <- 1:dimPred[match(dimDiff[d], attr(pred$Data, "dimensions"))]
      }
      indDimDiff<-as.matrix(expand.grid(indDimDiff))
      dimPermI <- match(attr(pred$Data, "dimensions")[setdiff(1:length(dimPred),match(dimDiff, attr(pred$Data, "dimensions")))],attr(obs$Data, "dimensions"))
      for (i in 1:dim(indDimDiff)[1]){
        if (is.null(window)){
          indTimePrd <- rep(list(bquote()), length(dimPred))
          indTimeSim <- rep(list(bquote()), length(dimFor))
          indTimePrd[match(dimDiff, attr(pred$Data, "dimensions"))] <- as.list(indDimDiff[i,])
          indTimeSim[match(dimDiff, attr(sim$Data, "dimensions"))] <- as.list(indDimDiff[i,])
          callPrd <- as.call(c(list(as.name("["),quote(pred$Data)), indTimePrd))
          callSim <- as.call(c(list(as.name("["),quote(sim$Data)), indTimeSim))
          #                              attrSim <- attr(sim$Data, "dimensions")
          if (!is.array(obs$Data)){
            F <- calibrateProj(obs$Data, eval(callPrd), eval(callSim), method = method, varcode = obs$Variable$varName, pr.threshold = threshold, scaling.type = scaling.type, extrapolate = extrapolation, theta=theta)
          }else{
            F <- calibrateProj(aperm(obs$Data, dimPermI), eval(callPrd), eval(callSim), method = method, varcode = obs$Variable$varName, pr.threshold = threshold, scaling.type = scaling.type, extrapolate = extrapolation, theta=theta)
          }
          indTimeSim <- rep(list(bquote()), length(dimFor))
          for (d in 1:length(dimFor)){
            indTimeSim[[d]] <- 1:dimFor[d]
          }
          indTimeSim[match(dimDiff, attr(sim$Data, "dimensions"))] <- as.list(indDimDiff[i,])
          indTimeSim<-as.matrix(expand.grid(indTimeSim))
          sim$Data[indTimeSim] <- F
          #                              attr(sim$Data, "dimensions") <- attrSim
        }else{
          for (j in 1:dim(dayList)[1]){
            indObs <- which(indDays == j)
            indObsWindow <- array(data = NaN, dim = c(window*length(indObs),1))
            for (d in 1:length(indObs)){
              indObsWindow[((d-1)*window+1):(d*window)] <- (indObs[d]-floor(window/2)):(indObs[d]+floor(window/2))
            }
            indObsWindow[which(indObsWindow <= 0)] <- 1
            indObsWindow[which(indObsWindow >  length(indDays))] <- length(indDays)
            indObsWindow <- unique(indObsWindow)
            indSim <- which(indDaysSim == j)
            indTimeObs <- rep(list(bquote()), length(dimObs))
            indTimePrd <- rep(list(bquote()), length(dimPred))
            indTimeSim <- rep(list(bquote()), length(dimFor))
            # Consideramos solo la ventana temporal
            indTimeObs[[obs.time.index]] <- indObsWindow
            indTimePrd[[pred.time.index]] <- indObsWindow
            indTimeSim[[pred.time.index]] <- indSim
            indTimePrd[match(dimDiff, attr(pred$Data, "dimensions"))] <- as.list(indDimDiff[i,])
            indTimeSim[match(dimDiff, attr(sim$Data, "dimensions"))] <- as.list(indDimDiff[i,])
            callObs <- as.call(c(list(as.name("["),quote(obs$Data)), indTimeObs))
            callPrd <- as.call(c(list(as.name("["),quote(pred$Data)), indTimePrd))
            callSim <- as.call(c(list(as.name("["),quote(sim$Data)), indTimeSim))
            if (!is.array(eval(callObs))){
              F <- calibrateProj(eval(callObs), eval(callPrd), eval(callSim), method = method, varcode = obs$Variable$varName, pr.threshold = threshold, scaling.type = scaling.type, extrapolate = extrapolation, theta=theta)
            }else{
              F <- calibrateProj(aperm(eval(callObs), dimPermI), eval(callPrd), eval(callSim), method = method, varcode = obs$Variable$varName, pr.threshold = threshold, scaling.type = scaling.type, extrapolate = extrapolation, theta=theta)
            }
            indTimeSim <- rep(list(bquote()), length(dimFor))
            for (d in 1:length(dimFor)){
              indTimeSim[[d]] <- 1:dimFor[d]
            }
            indTimeSim[[pred.time.index]] <- indSim
            indTimeSim[match(dimDiff, attr(sim$Data, "dimensions"))] <- as.list(indDimDiff[i,])
            indTimeSim<-as.matrix(expand.grid(indTimeSim))
            sim$Data[indTimeSim] <- F
          }
          #                              attr(sim$Data, "dimensions") <- attrSim
        }
      }
    }else{
      if (is.null(window)){
        obs.time.index <- grep("^time$", attr(obs$Data, "dimensions"))
        pred.time.index <- grep("^time$", attr(pred$Data, "dimensions"))
        pred.member.index <- grep("^member$", attr(pred$Data, "dimensions"))
        auxPerm <- c(pred.time.index, pred.member.index,setdiff(1:length(dimPred),c(pred.time.index, pred.member.index)))
        dimPermI <- match(attr(pred$Data, "dimensions")[setdiff(1:length(dimPred),match(dimDiff, attr(pred$Data, "dimensions")))],attr(obs$Data, "dimensions"))
        auxPrd <- array(data = rep(aperm(pred$Data, auxPerm),1), dim = c(dimPred[pred.member.index]*dimPred[pred.time.index],dimPred[setdiff(1:length(dimPred),c(pred.time.index, pred.member.index))]))
        auxSim <- array(data = rep(aperm(sim$Data, auxPerm),1), dim = c(dimFor[pred.member.index]*dimFor[pred.time.index],dimFor[setdiff(1:length(dimFor),c(pred.time.index, pred.member.index))]))
        auxObs <- array(data = NaN, dim = c(dimPred[pred.member.index]*dimPred[pred.time.index],dimPred[setdiff(1:length(dimPred),c(pred.time.index, pred.member.index))]))
        for (i in 1:dimPred[pred.member.index]){
          indMember <- ((i-1)*dimPred[pred.time.index]+1):(i*dimPred[pred.time.index])
          if (!is.array(obs$Data)){
            auxObs[indMember,,] <- obs$Data
          }else{
            auxObs[indMember,,] <- aperm(obs$Data, dimPermI)
          }
        }
        #                       attrSim <- attr(sim$Data, "dimensions")
        F <- calibrateProj(auxObs, auxPrd, auxSim, method = method, varcode = obs$Variable$varName, pr.threshold = threshold, scaling.type = scaling.type, extrapolate = extrapolation, theta=theta)
        sim$Data <- aperm(array(data = rep(F,1), dim = c(dimFor[pred.time.index],dimFor[pred.member.index],dimFor[setdiff(1:length(dimFor),c(pred.time.index, pred.member.index))])), auxPerm)
        #                        attr(sim$Data, "dimensions") <- attrSim
      }else{
        obs.time.index <- grep("^time$", attr(obs$Data, "dimensions"))
        pred.time.index <- grep("^time$", attr(pred$Data, "dimensions"))
        pred.member.index <- grep("^member$", attr(pred$Data, "dimensions"))
        auxPerm <- c(pred.time.index, pred.member.index,setdiff(1:length(dimPred),c(pred.time.index, pred.member.index)))
        dimPermI <- match(attr(pred$Data, "dimensions")[setdiff(1:length(dimPred),match(dimDiff, attr(pred$Data, "dimensions")))],attr(obs$Data, "dimensions"))
        for (j in 1:dim(dayList)[1]){
          indObs <- which(indDays == j)
          indObsWindow <- array(data = NaN, dim = c(window*length(indObs),1))
          for (d in 1:length(indObs)){
            indObsWindow[((d-1)*window+1):(d*window)] <- (indObs[d]-floor(window/2)):(indObs[d]+floor(window/2))
          }
          indObsWindow[which(indObsWindow <= 0)] <- 1
          indObsWindow[which(indObsWindow >  length(indDays))] <- length(indDays)
          indObsWindow <- unique(indObsWindow)
          indSim <- which(indDaysSim == j)
          indTimeObs <- rep(list(bquote()), length(dimObs))
          indTimePrd <- rep(list(bquote()), length(dimPred))
          indTimeSim <- rep(list(bquote()), length(dimFor))
          # Consideramos solo la ventana temporal
          indTimeObs[[obs.time.index]] <- indObsWindow
          indTimePrd[[pred.time.index]] <- indObsWindow
          indTimeSim[[pred.time.index]] <- indSim
          callObs <- as.call(c(list(as.name("["),quote(obs$Data)), indTimeObs))
          callPrd <- as.call(c(list(as.name("["),quote(pred$Data)), indTimePrd))
          callSim <- as.call(c(list(as.name("["),quote(sim$Data)), indTimeSim))
          auxPrd <- array(data = rep(aperm(eval(callPrd), auxPerm),1), dim = c(dimPred[pred.member.index]*length(indObsWindow),dimPred[setdiff(1:length(dimPred),c(pred.time.index, pred.member.index))]))
          auxSim <- array(data = rep(aperm(eval(callSim), auxPerm),1), dim = c(dimFor[pred.member.index]*length(indSim),dimFor[setdiff(1:length(dimFor),c(pred.time.index, pred.member.index))]))
          auxObs <- array(data = NaN, dim = c(dimPred[pred.member.index]*length(indObsWindow),dimPred[setdiff(1:length(dimPred),c(pred.time.index, pred.member.index))]))
          for (i in 1:dimPred[pred.member.index]){
            indMember <- ((i-1)*length(indObsWindow)+1):(i*length(indObsWindow))
            if (!is.array(eval(callObs))){
              auxObs[indMember,,] <- eval(callObs)
            }else{
              auxObs[indMember,,] <- aperm(eval(callObs), dimPermI)
            }
          }
          #                             attrSim <- attr(sim$Data, "dimensions")
          F <- calibrateProj(auxObs, auxPrd, auxSim, method = method, varcode = obs$Variable$varName, pr.threshold = threshold, scaling.type = scaling.type, extrapolate = extrapolation, theta=theta)
          F <- aperm(array(data = rep(F,1), dim = c(length(indSim),dimFor[pred.member.index],dimFor[setdiff(1:length(dimFor),c(pred.time.index, pred.member.index))])), auxPerm)
          indTimeSim <- rep(list(bquote()), length(dimFor))
          for (d in 1:length(dimFor)){
            indTimeSim[[d]] <- 1:dimFor[d]
          }
          indTimeSim[[pred.time.index]] <- indSim
          indTimeSim<-as.matrix(expand.grid(indTimeSim))
          sim$Data[indTimeSim] <- F
          #                              attr(sim$Data, "dimensions") <- attrSim
        }
      }
    }
  }
  if (any(grepl(obs$Variable$varName,c("pr","tp","precipitation","precip")))){
    attr(sim$Data, "threshold") <-  threshold
  }
  attr(sim$Data, "dimensions") <- attrSim
  attr(sim$Data, "correction") <-  method
  return(sim)
}
# End
################################################################################


#' @title Bias correction methods
#' @description Implementation of several standard bias correction methods
#'
#' @template templateObsPredSim
#' @param method method applied. Current accepted values are \code{"eqm"}, \code{"delta"},
#' \code{"scaling"}, \code{"gqm"} and \code{"gpqm"}. See details.
#' @param varcode Variable code. This is not the variable itself to be corrected, but
#' rather it referes to its nature and distributional properties. For instance, \code{"tas"} applies for
#' temperature-like variables (i.e.: unbounded gaussian variables that can take both negative and positive values),
#' \code{"hurs"} for relative-humidity-like variables (i.e., bounded by both sides, with a maximum value of 1 -100-
#' and a minimum of zero, with no negatives), \code{"wss"} for wind-like variables (no possible values below zero,
#' but without an upper bound) and \code{"pr"}, specifically for precipitation.
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#' \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm).
#' @param scaling.type Character indicating the type of the scaling method. Options are \code{"additive"} (default)
#' or \code{"multiplicative"}. This argument is ignored if \code{"scaling"} is not 
#' selected as the bias correction method.
#' @param extrapolate Character indicating the extrapolation method to be applied to correct values in  
#' \code{"sim"} that are out of the range of \code{"pred"}. Extrapolation is applied only to the \code{"eqm"} method, 
#' thus, this argument is ignored if other bias correction method is selected. Default is \code{"no"} (do not extrapolate).
#' @param theta numeric indicating  upper threshold (and lower for the left tail of the distributions, if needed) 
#' above which precipitation (temperature) values are fitted to a Generalized Pareto Distribution (GPD). 
#' Values below this threshold are fitted to a gamma (normal) distribution. By default, 'theta' is the 95th 
#' percentile (5th percentile for the left tail). Only for \code{"gpqm"} method.
#' 
#' 
#@author S. Herrera \email{sixto@@predictia.es}
#' @export
#' @family downscaling
#' @family calibration
#' @importFrom MASS fitdistr
#' @importFrom evd fpot
#' @importFrom evd qgpd
#' @importFrom evd pgpd
#' @keywords internal
#'


calibrateProj <- function (obs, pred, sim, method = c("eqm", "delta", "scaling", "gqm", "gpqm"), varcode = c("tas", "hurs", "tp", "pr", "wss"), pr.threshold = 1, scaling.type = scaling.type, extrapolate = c("no", "constant"), theta = .95) {
  d1 <- FALSE
  d2 <- FALSE
  if (!is.array(obs)){
    obs <- array(data = obs, dim = c(length(obs),1,1))
    pred <- array(data = pred, dim = c(length(pred),1,1))
    sim <- array(data = sim, dim = c(length(sim),1,1))
    d1 <- TRUE
  }else{
    if (length(dim(obs))==2){
      obs <- array(data = obs, dim = c(dim(obs),dim(obs)[2]))
      pred <- array(data = pred, dim = c(dim(pred),dim(pred)[2]))
      sim <- array(data = sim, dim = c(dim(sim),dim(sim)[2]))
      d2 <- TRUE
    }
  }
  if (any(grepl(varcode,c("pr","tp","precipitation","precip")))) {
    threshold<-pr.threshold
    nP<-matrix(data = NA, ncol=dim(pred)[3], nrow=dim(pred)[2])
    Pth<-matrix(data = NA, ncol=dim(pred)[3], nrow=dim(pred)[2])
    for (i in 1:dim(pred)[2]){
      for (j in 1:dim(pred)[3]){
        nP[i,j]<-sum(as.double(obs[,i,j]<=threshold & !is.na(obs[,i,j])), na.rm = TRUE)
        if (nP[i,j]>=0 & nP[i,j]<dim(obs)[1]){
          ix<-sort(pred[,i,j], decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
          Ps<-sort(pred[,i,j], decreasing = FALSE, na.last = NA)
          Pth[i,j]<-Ps[nP[i,j]+1]
          if (Ps[nP[i,j]+1]<=threshold){
            Os<-sort(obs[,i,j], decreasing = FALSE, na.last = NA)
            ind<-which(Ps > threshold & !is.na(Ps))
            if (length(ind)==0){
              ind <- max(which(!is.na(Ps)))
              ind <- min(c(length(Os),ind))
            }else{
              ind<-min(which(Ps > threshold & !is.na(Ps)))
            }
            # [Shape parameter Scale parameter]
            if (length(unique(Os[(nP[i,j]+1):ind]))<6){
              Ps[(nP[i,j]+1):ind] <- mean(Os[(nP[i,j]+1):ind], na.rm = TRUE)
            }else{
              
              auxOs <- Os[(nP[i,j]+1):ind]
              auxOs <- auxOs[which(!is.na(auxOs))]
              auxGamma <- fitdistr(auxOs,"gamma")
              
              Ps[(nP[i,j]+1):ind]<-rgamma(ind-nP[i,j], auxGamma$estimate[1], rate = auxGamma$estimate[2])
            }
            Ps<-sort(Ps, decreasing = FALSE, na.last = NA)
          }
          if (nP[i,j]>0){
            ind<-min(nP[i,j],dim(pred)[1])
            Ps[1:ind]<-0
          }
          pred[ix,i,j]<-Ps
        }else{
          if (nP[i,j]==dim(obs)[1]){
            ix<-sort(pred[,i,j], decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
            Ps<-sort(pred[,i,j], decreasing = FALSE, na.last = NA)
            Pth[i,j]<-Ps[nP[i,j]]
            ind<-min(nP[i,j],dim(pred)[1])
            Ps[1:ind]<-0
            pred[ix,i,j]<-Ps
          }
        }
      }
    }
  }
  if (method == "delta") {
    if (dim(sim)[1]!=dim(obs)[1]){
      stop("sim and obs should have the same dimensions")
    }else{
      simMean <- apply(sim, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
      prdMean <- apply(pred, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
      for (k in 1:min(c(dim(sim)[1],dim(obs)[1]))) {
        sim[k,,]<-obs[k,,]-prdMean+simMean
      }
    }
  }
  if (method == "scaling") {
    if (scaling.type == "additive"){
      obsMean <- apply(obs, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
      prdMean <- apply(pred, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
      for (k in 1:dim(sim)[1]) {
        sim[k,,]<-sim[k,,]-prdMean+obsMean
      }
    }else if(scaling.type == "multiplicative"){
      obsMean <- apply(obs, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
      prdMean <- apply(pred, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
      for (k in 1:dim(sim)[1]) {
        sim[k,,]<-(sim[k,,]/prdMean)*obsMean
      }
    }
  }
  
  if (method == "eqm") {
    if (!any(grepl(varcode,c("pr","tp","precipitation","precip")))){
      for (i in 1:dim(sim)[2]) {
        for (j in 1:dim(sim)[3]) {
          if (any(!is.na(pred[,i,j])) & any(!is.na(obs[,i,j]))) {
            ePrd <- ecdf(pred[,i,j])
            ################################################################################################3           
            if(extrapolate == "constant"){
              exupsim <- which(sim[,i,j] > max(range(pred[,i,j], na.rm = TRUE))) 
              exdwnsim <- which(sim[,i,j] < min(range(pred[,i,j], na.rm = TRUE)))  
              ex <- c(exupsim,exdwnsim)
              exC <- setdiff(1:dim(sim)[1],ex)
              if (length(exupsim)>0){
                extmin <- max(pred[,i,j], na.rm = TRUE)
                dif <- extmin - max(obs[,i,j], na.rm = TRUE)
                sim[exupsim,i,j] <- sim[exupsim,i,j] - dif
              }
              if (length(exdwnsim)>0){
                extmin <- min(pred[,i,j], na.rm = TRUE)
                dif <- extmin - min(obs[,i,j], na.rm = TRUE)
                sim[exdwnsim,i,j] <- sim[exdwnsim,i,j] - dif
              }
              sim[exC,i,j] <- quantile(obs[,i,j], probs = ePrd(sim[exC,i,j]), na.rm = TRUE, type = 4)
            }else{
              sim[,i,j] <- quantile(obs[,i,j], probs = ePrd(sim[,i,j]), na.rm = TRUE, type = 4)
            }
            
            ################################################################################################3
          }
        }
      }
    } else {
      for (i in 1:dim(sim)[2]) {
        for (j in 1:dim(sim)[3]) {
          if (any(!is.na(pred[,i,j])) & any(!is.na(obs[,i,j]))){
            if (length(which(pred[,i,j]>Pth[i,j]))>0){ 
              ##
              ePrd<-ecdf(pred[which(pred[,i,j]>Pth[i,j]),i,j]) 
              ##
              noRain<-which(sim[,i,j]<=Pth[i,j] & !is.na(sim[,i,j]))
              rain<-which(sim[,i,j]>Pth[i,j] & !is.na(sim[,i,j]))
              drizzle<-which(sim[,i,j]>Pth[i,j] & sim[,i,j] <= min(pred[which(pred[,i,j]>Pth[i,j]),i,j], na.rm = TRUE) & !is.na(sim[,i,j]))
              if (length(rain)>0){
                eFrc<-ecdf(sim[rain,i,j]) #eFrc??
                ##############################################################################3
                if(extrapolate == "constant"){
                  exupsim <- which(sim[rain,i,j] > max(range(pred[,i,j], na.rm = TRUE))) #
                  exdwnsim <- which(sim[rain,i,j] < min(range(pred[,i,j], na.rm = TRUE)))  
                  ex <- c(exupsim,exdwnsim)
                  exC <- setdiff(1:length(rain),ex)
                  if (length(exupsim)>0){
                    extmin <- max(pred[,i,j], na.rm = TRUE)
                    dif <- extmin - max(obs[which(obs[,i,j] > threshold & !is.na(obs[,i,j])),i,j], na.rm = TRUE)
                    sim[rain[exupsim],i,j] <- sim[rain[exupsim],i,j] - dif
                  }
                  if (length(exdwnsim)>0){
                    extmin <- min(pred[,i,j], na.rm = TRUE)
                    dif <- extmin - min(obs[which(obs[,i,j] > threshold & !is.na(obs[,i,j])),i,j], na.rm = TRUE)
                    sim[rain[exdwnsim],i,j] <- sim[rain[exdwnsim],i,j] - dif
                  }
                  sim[rain[exC],i,j] <- quantile(obs[which(obs[,i,j] > threshold & !is.na(obs[,i,j])),i,j], probs = ePrd(sim[rain[exC],i,j]), na.rm = TRUE, type = 4)
                }else{
                  sim[rain,i,j]<-quantile(obs[which(obs[,i,j]>threshold & !is.na(obs[,i,j])),i,j], probs = ePrd(sim[rain,i,j]), na.rm = TRUE, type = 4)
                }
                ################################################################################################3
              }
              if (length(drizzle)>0){
                sim[drizzle,i,j]<-quantile(sim[which(sim[,i,j]>min(pred[which(pred[,i,j]>Pth[i,j]),i,j], na.rm = TRUE) & !is.na(sim[,i,j])),i,j], probs = eFrc(sim[drizzle,i,j]), na.rm = TRUE, type = 4)
              }
              sim[noRain,i,j]<-0
            }else{
              noRain<-which(sim[,i,j]<=Pth[i,j] & !is.na(sim[,i,j]))
              rain<-which(sim[,i,j]>Pth[i,j] & !is.na(sim[,i,j]))
              if (length(rain)>0){
                eFrc<-ecdf(sim[rain,i,j])
                sim[rain,i,j]<-quantile(obs[which(obs[,i,j] > threshold & !is.na(obs[,i,j])),i,j], probs = eFrc(sim[rain,i,j]), na.rm = TRUE, type = 4)
              }
              sim[noRain,i,j]<-0
            }
          }
        }
      }
    }
  }
  if (method == "gqm" & any(grepl(varcode,c("pr","tp","precipitation","precip")))) {
    for (i in 1:dim(sim)[2]){
      for (j in 1:dim(sim)[3]){
#        if (nP[i,j]>0){
        if (nP[i,j]<dim(obs)[1]){
          ind<-which(obs[,i,j]>threshold & !is.na(obs[,i,j]))
          obsGamma<-fitdistr(obs[ind,i,j],"gamma")
          ind<-which(pred[,i,j]>0 & !is.na(pred[,i,j]))
          prdGamma<-fitdistr(pred[ind,i,j],"gamma")
          rain<-which(sim[,i,j]>Pth[i,j] & !is.na(sim[,i,j]))
          noRain<-which(sim[,i,j]<=Pth[i,j] & !is.na(sim[,i,j]))
          auxF<-pgamma(sim[rain,i,j],prdGamma$estimate[1], rate = prdGamma$estimate[2])
          sim[rain,i,j]<-qgamma(auxF,obsGamma$estimate[1], rate = obsGamma$estimate[2])
          sim[noRain,i,j]<-0
          warningNoRain <- FALSE
        }else{
          warningNoRain <- TRUE
        }
      }
    }
    if (warningNoRain){
      warning("There is at least one location without rainfall above the threshold.\n In this (these) location(s) none bias correction has been applied.")
    }
  }
  if (method == "gpqm" & any(grepl(varcode,c("pr","tp","precipitation","precip")))) {
    for (i in 1:dim(sim)[2]){
      for (j in 1:dim(sim)[3]){
        #        if (nP[i,j]>0){
        if (nP[i,j]<dim(obs)[1]){
          ind<-which(obs[,i,j]>threshold & !is.na(obs[,i,j]))
          indgamma <- ind[which(obs[ind,i,j] < quantile(obs[ind,i,j], theta))]
          indpareto <- ind[which(obs[ind,i,j] >= quantile(obs[ind,i,j], theta))]
          obsGQM <- fitdistr(obs[indgamma,i,j],"gamma")
          obsGQM2 <- fpot(obs[indpareto,i,j], quantile(obs[ind,i,j], theta), "gpd", std.err = FALSE)
          ind <- which(pred[,i,j] > 0 & !is.na(pred[,i,j]))
          indgamma <- ind[which(pred[ind,i,j]<quantile(pred[ind,i,j],theta))]
          indpareto <-ind[which(pred[ind,i,j]>=quantile(pred[ind,i,j], theta))]
          prdGQM <- fitdistr(pred[indgamma,i,j], "gamma")
          prdGQM2 <- fpot(pred[indpareto,i,j], quantile(pred[ind,i,j], theta), "gpd", std.err = FALSE)
          rain <- which(sim[,i,j] > Pth[i,j] & !is.na(sim[,i,j]))
          noRain<-which(sim[,i,j] <= Pth[i,j] & !is.na(sim[,i,j]))
          indgammasim <- rain[which(sim[rain,i,j] < quantile(pred[ind,i,j], theta))]
          indparetosim <- rain[which(sim[rain,i,j] >= quantile(pred[ind,i,j], theta))]
          auxF <- pgamma(sim[indgammasim,i,j],prdGQM$estimate[1], rate = prdGQM$estimate[2])
          auxF2 <- pgpd(sim[indparetosim,i,j],loc = 0, scale = prdGQM2$estimate[1], shape = prdGQM2$estimate[2])
          sim[indgammasim,i,j] <- qgamma(auxF, obsGQM$estimate[1], rate = obsGQM$estimate[2])
          sim[indparetosim[which(auxF2<1)],i,j] <- qgpd(auxF2[which(auxF2 < 1)], loc = 0, scale = obsGQM2$estimate[1], shape = obsGQM2$estimate[2])
          sim[indparetosim[which(auxF2==1)],i,j] <- max(obs[indpareto,i,j], na.rm = TRUE)
          sim[noRain,i,j]<-0
          warningNoRain <- FALSE
        }else{
          warningNoRain <- TRUE
        }
      }
    }
    if (warningNoRain){
      warning("There is at least one location without rainfall above the threshold.\n In this (these) location(s) none bias correction has been applied.")
    }
  }
  if (d1){
    sim <- sim[,1,1]
  }else{
    if (d2){
      sim <- sim[,,1]
    }
  }
  return(sim)
}

# End
