#' @title Bias correction methods
#' @description Implementation of several standard bias correction methods
#' @author S. Herrera \email{sixto@@predictia.es}
#' @export
#' 

biasCorrection <- function (obs, pred, sim, method = c("qqmap", "delta", "scaling", "unbiasing", "piani"), pr.threshold = 1) {
      
      method<-match.arg(method, choices = c("qqmap", "delta", "unbiasing", "piani", "scaling"))
      threshold<-pr.threshold
      dimObs<-dim(obs$Data)
      dimPred<-dim(pred$Data)
      dimFor<-dim(sim$Data)
      dimDiff<-NULL
      dimPerm <- 1:length(attr(pred$Data, "dimensions"))
      if (length(setdiff(attr(pred$Data, "dimensions"),attr(obs$Data, "dimensions")))>0){
            dimDiff<-setdiff(attr(pred$Data, "dimensions"),attr(obs$Data, "dimensions"))
            for (k in 1:length(dimDiff)) {
                  dimPerm[which(grepl(dimDiff[k],attr(pred$Data, "dimensions")))]<-length(attr(obs$Data, "dimensions"))+k
            }
      }
      dimPerm[which(dimPerm<=length(dim(obs$Data)))]<-  match(attr(obs$Data, "dimensions"), attr(pred$Data, "dimensions"))
      
      if (length(dimDiff)==0){
            F <- calibrateProj(obs$Data, aperm(pred$Data,dimPerm), aperm(sim$Data,dimPerm), method = method, varcode = obs$Variable$varName, pr.threshold = threshold)
            sim$Data<-aperm(F,c(dimPerm[2:length(dim(obs$Data))],dimPerm[1]))
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
      #       if (length(dimDiff)==1){
      #             for (n in 1:dim(sim$Data)[length(dim(obs$Data))+1]){
      #                   F1 <- calibrateProj(obs$Data, P[,,,n], F[,,,n], method = method, varcode = obs$Variable$varName, pr.threshold = threshold)
      #                   #    sim$Data[,,,n]<-aperm(F1,c(dimPerm[2:length(dim(obs$Data))],dimPerm[1],dimPerm[(length(dim(obs$Data))+1):length(dim(pred$Data))]))
      #                   sim$Data[,,,n]<-aperm(F1,c(dimPerm[2:length(dim(obs$Data))],dimPerm[1]))
      #             }
      #       }
      #       if (length(dimDiff)==2){
      #             for (n in 1:dim(sim$Data)[length(dim(obs$Data))+1]){
      #                   for (m in 1:dim(sim$Data)[length(dim(obs$Data))+2]){
      #                         F1 <- calibrateProj(obs$Data, P[,,,n,m], F[,,,n,m], method = method, varcode = obs$Variable$varName, pr.threshold = threshold)
      #                         sim$Data[,,,n,m]<-aperm(F1,c(dimPerm[2:length(dim(obs$Data))],dimPerm[1]))
      #                   }
      #             }
      #       }
      #       if (length(dimDiff)==3){
      #             for (n in 1:dim(sim$Data)[length(dim(obs$Data))+1]){
      #                   for (m in 1:dim(sim$Data)[length(dim(obs$Data))+2]){
      #                         for (t in 1:dim(sim$Data)[length(dim(obs$Data))+3]){
      #                               F1 <- calibrateProj(obs$Data, P[,,,n,m,t], F[,,,n,m,t], method = method, varcode = obs$Variable$varName, pr.threshold = threshold)
      #                               sim$Data[,,,n,m,t]<-aperm(F1,c(dimPerm[2:length(dim(obs$Data))],dimPerm[1]))
      #                         }
      #                   }
      #             }
      #       }
      if (any(grepl(obs$Variable$varName,c("pr","tp","precipitation","precip")))){
            attr(sim$Data, "threshold") <-  threshold
      }
      attr(sim$Data, "correction") <-  method
      return(sim)
}
