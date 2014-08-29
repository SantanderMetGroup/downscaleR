#' @title Bias correction methods
#' @description Implementation of several standard bias correction methods
#' 
#' @template templateObsPredSim
#' @param method method applied. Current accepted values are \code{"qqmap"}, \code{"delta"},
#'  \code{"scaling"}, \code{"unbiasing"} and \code{"piani"}. See details.
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm).
#'  
#'  @details ~~ Details on the different methods here
#'  
#'  @return A calibrated object of the same spatio-temporal extent of the input field

#' @author S. Herrera \email{sixto@@predictia.es}
#' @export
#' @family downscaling
#' @family calibration
#' @references ~References to the different methods
#' 

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
