#' @title Bias correction methods
#' @description Implementation of several standard bias correction methods
#' 
#' @template templateObsPredSim
#' @param method method applied. Current accepted values are \code{"qqmap"}, \code{"delta"},
#'  \code{"scaling"}, \code{"unbiasing"} and \code{"piani"}. See details.
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm).
#'  
#' @details
#' 
#' The methods available are \code{"qqmap"}, \code{"delta"}, \code{"unbiasing"}, 
#' \code{"scaling"} and \code{"Piani"} (the latter used only for precipitation).
#' Next we make a brief description of each method:
#' 
#' \strong{Delta}
#' 
#' This method consists on adding to the observations the mean change signal (delta method).
#' This method is applicable to any kind of variable but it is preferable to avoid it for bounded variables
#'  (e.g. precipitation, wind speed, etc.) because values out of the variable range could be obtained (e.g. negative wind speeds...).
#' 
#' \strong{Unbiasing}
#' 
#' This correction consists on adding to the simulation the mean diference between the observations 
#' and the simulation in the train period. This method is preferably applicable to unbounded variables (e.g. temperature).
#' 
#' \strong{Scaling}
#' 
#' This method consists on scaling the simulation by the quotient of the mean observed and simulated in the train period. 
#' This method is preferably applicable to variables with a lower bound, such as precipitation, because it preserves the 
#' frequency also.
#' 
#' \strong{Quantile-Quantile Mapping (qqmap)}
#' 
#' This is a very extended bias correction method which consists on calibrating the simulated Cumulative Distribution Function (CDF) 
#' by adding to the observed quantiles both the mean delta change and the individual delta changes in the corresponding quantiles. 
#' This method is applicable to any kind of variable.
#' 
#' \strong{Piani}
#' 
#' This method is described in Piani et al. 2010 and is applicable only to precipitation. It is based on the initial assumption that both observed
#'  and simulated intensity distributions are well approximated by the gamma distribution, therefore is a parametric q-q map 
#'  that uses the theorical instead of the empirical distribution. 
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
#' }
#' @author S. Herrera \email{sixto@@predictia.es}
#' @export
#' @examples \dontrun{
#' # These are the paths to the package built-in GSN and NCEP datasets 
#' gsn.data.dir <- file.path(find.package("downscaleR"), "datasets/observations/GSN_Iberia")
#' ncep.data.dir <- file.path(find.package("downscaleR"), "datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.ncml")
#' gsn.inv <- dataInventory(gsn.data.dir)
#' ncep.inv <- dataInventory(ncep.data.dir)
#' str(gsn.inv)
#' str(ncep.inv)
#' # Load precipitation for boreal winter (DJF) in the train (1991-2000) and test (2001-2010) periods, for the observations (GSN_Iberia) and the Iberia_NCEP datasets
#' obs <- loadStationData(dataset = gsn.data.dir, var="precip", lonLim = c(-12,10), latLim = c(33,47), season=c(12,1,2), years = 1991:2000)
#' prd <- loadGridData(ncep.data.dir, var = "tp", lonLim = c(-12,10), latLim = c(33,47), season = c(12,1,2), years = 1991:2000)
#' sim <- loadGridData(ncep.data.dir, var = "tp", lonLim = c(-12,10), latLim = c(33,47), season = c(12,1,2), years = 2001:2010)
#' # Interpolation of the observations onto the grid of model: we use the method "nearest" and the getGrid function to ensure spatial consistency:
#' obs <- interpGridData(obs, new.grid = getGrid(prd), method = "nearest")
#' # Apply the bias correction method:
#' simBC <- biasCorrection (obs, prd, sim, method = "qqmap", pr.threshold = 1) # qq-mapping
#' par(mfrow = c(1,2))
#' plotMeanField(sim)
#' plotMeanField(simBC)
#' }


biasCorrection <- function (obs, pred, sim, method = c("qqmap", "delta", "scaling", "unbiasing", "piani"), pr.threshold = 1) {
      method<-match.arg(method, choices = c("qqmap", "delta", "unbiasing", "piani", "scaling"))
      threshold<-pr.threshold
      dimObs<-dim(obs$Data)
      dimPred<-dim(pred$Data)
      dimFor<-dim(sim$Data)
      dimDiff<-NULL
      dimPerm <- 1:length(attr(pred$Data, "dimensions"))
      dimPerm[which(dimPerm<=length(dim(obs$Data)))]<- match(attr(obs$Data, "dimensions"), attr(pred$Data, "dimensions"))
      if (length(setdiff(attr(pred$Data, "dimensions"),attr(obs$Data, "dimensions")))>0){
            dimDiff<-setdiff(attr(pred$Data, "dimensions"),attr(obs$Data, "dimensions"))
            for (k in 1:length(dimDiff)) {
                  indDiff <- which(grepl(dimDiff[k],attr(pred$Data, "dimensions")))
                  dimPerm[length(attr(obs$Data, "dimensions"))+k] <- indDiff
            }
      }
      if (length(dimDiff)==0){
            F <- calibrateProj(obs$Data, aperm(pred$Data,dimPerm), aperm(sim$Data,dimPerm), method = method, varcode = obs$Variable$varName, pr.threshold = threshold)
            if (!any(dimPerm != 1:length(attr(pred$Data, "dimensions")))){
                  sim$Data<-F
            }else{
                  sim$Data<-aperm(F,c(dimPerm[2:length(dim(obs$Data))],dimPerm[1]))
            }
      }else{
            indDimDiff <- rep(list(bquote()), length(dimDiff))
            for (d in 1:length(dimDiff)){
                  indDimDiff[[d]] <- 1:dimPred[match(dimDiff[d], attr(pred$Data, "dimensions"))]
            }
            indDimDiff<-as.matrix(expand.grid(indDimDiff))
            dimPermI <- match(attr(pred$Data, "dimensions")[setdiff(1:length(dimPred),match(dimDiff, attr(pred$Data, "dimensions")))],attr(obs$Data, "dimensions"))
            for (i in 1:dim(indDimDiff)[1]){
                  indTimePrd <- rep(list(bquote()), length(dimPred))
                  indTimeSim <- rep(list(bquote()), length(dimFor))
                  indTimePrd[match(dimDiff, attr(pred$Data, "dimensions"))] <- as.list(indDimDiff[i,])
                  indTimeSim[match(dimDiff, attr(sim$Data, "dimensions"))] <- as.list(indDimDiff[i,])
                  callPrd <- as.call(c(list(as.name("["),quote(pred$Data)), indTimePrd))
                  callSim <- as.call(c(list(as.name("["),quote(sim$Data)), indTimeSim))
                  F <- calibrateProj(aperm(obs$Data, dimPermI), eval(callPrd), eval(callSim), method = method, varcode = obs$Variable$varName, pr.threshold = threshold)
                  indTimeSim <- rep(list(bquote()), length(dimFor))
                  for (d in 1:length(dimFor)){
                        indTimeSim[[d]] <- 1:dimFor[d]
                  }
                  indTimeSim[match(dimDiff, attr(sim$Data, "dimensions"))] <- as.list(indDimDiff[i,])
                  indTimeSim<-as.matrix(expand.grid(indTimeSim))
                  sim$Data[indTimeSim] <- F
            }
            
      }
      if (any(grepl(obs$Variable$varName,c("pr","tp","precipitation","precip")))){
            attr(sim$Data, "threshold") <-  threshold
      }
      attr(sim$Data, "correction") <-  method
      return(sim)
}
# End
################################################################################


#' @title Bias correction methods
#' @description Implementation of several standard bias correction methods
#'
#' @template templateObsPredSim
#' @param method method applied. Current accepted values are \code{"qqmap"}, \code{"delta"},
#' \code{"scaling"}, \code{"unbiasing"} and \code{"piani"}. See details.
#' @param varcode Variable code. This is not the variable itself to be corrected, but
#' rather it referes to its nature and distributional properties. For instance, \code{"tas"} applies for
#' temperature-like variables (i.e.: unbounded gaussian variables that can take both negative and positive values),
#' \code{"hurs"} for relative-humidity-like variables (i.e., bounded by both sides, with a maximum value of 1 -100-
#' and a minimum of zero, with no negatives), \code{"wss"} for wind-like variables (no possible values below zero,
#' but without an upper bound) and \code{"pr"}, specifically for precipitation.
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#' \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm).
#' @author S. Herrera \email{sixto@@predictia.es}
#' @export
#' @family downscaling
#' @family calibration
#' @importFrom MASS fitdistr
#' @keywords internal
#'
calibrateProj <- function (obs, pred, sim, method = c("qqmap", "delta", "scaling", "unbiasing", "piani"), varcode = c("tas", "hurs", "pr", "wss"), pr.threshold = 1) {
      if (varcode == "pr") {
            threshold<-pr.threshold
            nP<-matrix(data = NA, ncol=dim(pred)[3], nrow=dim(pred)[2])
            Pth<-matrix(data = NA, ncol=dim(pred)[3], nrow=dim(pred)[2])
            for (i in 1:dim(pred)[2]){
                  for (j in 1:dim(pred)[3]){
                        nP[i,j]<-sum(as.double(obs[,i,j]<=threshold & !is.na(obs[,i,j])), na.rm = TRUE)
                        if (nP[i,j]>0){
                              ix<-sort(pred[,i,j], decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
                              Ps<-sort(pred[,i,j], decreasing = FALSE, na.last = NA)
                              Pth[i,j]<-Ps[nP[i,j]+1]
                              if (Ps[nP[i,j]+1]<=threshold){
                                    ind<-min(which(Ps > threshold & !is.na(Ps)))
                                    # [Shape parameter Scale parameter]
                                    Os<-sort(obs[,i,j], decreasing = FALSE, na.last = NA)
                                    auxGamma<-fitdistr(Os[(nP[i,j]+1):ind],"gamma")
                                    Ps[(nP[i,j]+1):ind]<-rgamma(ind-nP[i,j], auxGamma$estimate[1], rate = auxGamma$estimate[2])
                                    Ps<-sort(Ps, decreasing = FALSE, na.last = NA)
                              }
                              ind<-min(nP[i,j],dim(pred)[1])
                              Ps[1:ind]<-0
                              pred[ix,i,j]<-Ps
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
      if (method == "unbiasing") {
            obsMean <- apply(obs, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
            prdMean <- apply(pred, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
            for (k in 1:dim(sim)[1]) {
                  sim[k,,]<-sim[k,,]-prdMean+obsMean
            }
      }
      if (method == "scaling") {
            obsMean <- apply(obs, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
            prdMean <- apply(pred, FUN = mean, MARGIN = c(2,3), na.rm = TRUE)
            for (k in 1:dim(sim)[1]) {
                  sim[k,,]<-(sim[k,,]/prdMean)*obsMean
            }
      }
      if (method == "qqmap") {
            if (varcode != "pr") {      
                  for (i in 1:dim(sim)[2]) {
                        for (j in 1:dim(sim)[3]) {
                              if (any(!is.na(pred[,i,j])) & any(!is.na(obs[,i,j]))) {
                                    ePrd<-ecdf(pred[,i,j])
                                    sim[,i,j]<-quantile(obs[,i,j], probs = ePrd(sim[,i,j]), na.rm = TRUE, type = 4)
                              }
                        }
                  }
            } else {
                  for (i in 1:dim(sim)[2]) {
                        for (j in 1:dim(sim)[3]) {
                              if (any(!is.na(pred[,i,j])) & any(!is.na(obs[,i,j]))){
                                    ePrd<-ecdf(pred[which(pred[,i,j]>Pth[i,j]),i,j])
                                    noRain<-which(sim[,i,j]<=Pth[i,j] & !is.na(sim[,i,j]))
                                    rain<-which(sim[,i,j]>Pth[i,j] & !is.na(sim[,i,j]))
                                    drizzle<-which(sim[,i,j]>Pth[i,j] & sim[,i,j]<=min(pred[which(pred[,i,j]>Pth[i,j]),i,j], na.rm = TRUE) & !is.na(sim[,i,j]))
                                    eFrc<-ecdf(sim[rain,i,j])
                                    sim[drizzle,i,j]<-quantile(sim[which(sim[,i,j]>min(pred[which(pred[,i,j]>Pth[i,j]),i,j], na.rm = TRUE) & !is.na(sim[,i,j])),i,j], probs = eFrc(sim[drizzle,i,j]), na.rm = TRUE, type = 4)
                                    sim[rain,i,j]<-quantile(obs[which(obs[,i,j]>threshold & !is.na(obs[,i,j])),i,j], probs = ePrd(sim[rain,i,j]), na.rm = TRUE, type = 4)
                                    sim[noRain,i,j]<-0
                              }
                        }
                  }
            }
      }
      if (method == "piani" & varcode == "pr") {
            for (i in 1:dim(sim)[2]){
                  for (j in 1:dim(sim)[3]){
                        if (nP[i,j]>0){
                              ind<-which(obs[,i,j]>threshold & !is.na(obs[,i,j]))
                              obsGamma<-fitdistr(obs[ind,i,j],"gamma")
                              ind<-which(pred[,i,j]>0 & !is.na(pred[,i,j]))
                              prdGamma<-fitdistr(pred[ind,i,j],"gamma")
                              rain<-which(sim[,i,j]>Pth[i,j] & !is.na(sim[,i,j]))
                              noRain<-which(sim[,i,j]<=Pth[i,j] & !is.na(sim[,i,j]))
                              auxF<-pgamma(sim[rain,i,j],prdGamma$estimate[1], rate = prdGamma$estimate[2])
                              sim[rain,i,j]<-qgamma(auxF,obsGamma$estimate[1], rate = obsGamma$estimate[2])
                              sim[noRain,i,j]<-0
                        }
                  }
            }
      }
      return(sim)
}
# End

