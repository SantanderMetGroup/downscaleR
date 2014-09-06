#' @title Bias correction ISI-MIP method
#' @description Implementation of the ISI-MIP methodology
#'
#' @template templateObsPredSim
#' @param pr.threshold The minimum value that is considered as a non-zero precipitation. Ignored for
#'  \code{varcode} values different from \code{"pr"}. Default to 1 (assuming mm).
#' @details
#' 
#' The methods available are qqmap, delta, unbiasing, scaling and Piani (only precipitation).
#' 
#' \strong{ISI-MIP}
#' 
#' Recently, Hempel et al.2013 proposed a new bias correction methodology within the ISI-MIP Project, 
#' the first Inter-Sectoral Impact Model Intercomparison Project, funded by the German Federal Ministry 
#' of Education and Research (BMBF). This method has been developed to preserve the change signal (trend, 
#' climate change signal, etc.) and can be applied to several variables (precipitation, mean, maximum and 
#' minimum temperature, windspeed and eastward/northward components, radiation, pressure and humidity). 
#' The main difference with the rest of bias correction methods included in the \code{\link{biasCorrection}} 
#' function is that the ISI-MIP method includes dependencies between some variables. That is, to correct some 
#' of the variables (maximum/minimum temperatures and eastward/northward wind components) others are needed 
#' (mean temperature and windspeed).
#' 
#' sim <- isimip(obs, pred, sim) # Temperature or other variable
#'
#' sim <- isimip(obs, pred, sim, pr.threshold = threshold) # In the case of precipitation we should include the threshold considered of wet/dry days
#'
#' @seealso \code{\link{biasCorrection}} for details on other standard methods for bias correction
#'  
#' @return A calibrated object of the same spatio-temporal extent of the input field
#' 
#' @family downscaling
#' 
#' @references
#' 
#' \itemize{
#' \item Hempel, S., Frieler, K., Warszawski, L., Schewe, J., and Piontek, F. (2013) A trend-preserving bias correction: the ISI-MIP approach, Earth Syst. Dynam., 4, 219-236
#' }
#' @author S. Herrera \email{sixto@@predictia.es}
#' @export
#' @examples \dontrun{
#' # These are the paths to the package built-in GSN and NCEP datasets 
#' gsn.data.dir <- file.path(find.package("downscaleR"), "datasets/observations/GSN_Iberia")
#' ncep.data.dir <- file.path(find.package("downscaleR"), "datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.ncml")
#' # Data inventories provides a quick overview of the available data
#' gsn.inv <- dataInventory(gsn.data.dir)
#' ncep.inv <- dataInventory(ncep.data.dir)
#' str(gsn.inv)
#' str(ncep.inv)
#' # Load precipitation for boreal winter (DJF) in the train (1991-2000) and test (2001-2010) periods,
#' # for the observations (GSN_Iberia) and the Iberia_NCEP datasets
#' obs <- loadStationData(dataset = gsn.data.dir, var="precip", lonLim = c(-12,10), latLim = c(33,47), season=c(12,1,2), years = 1991:2000)
#' prd <- loadGridData(ncep.data.dir, var = "tp", lonLim = c(-12,10), latLim = c(33,47), season = c(12,1,2), years = 1991:2000)
#' sim <- loadGridData(ncep.data.dir, var = "tp", lonLim = c(-12,10), latLim = c(33,47), season = c(12,1,2), years = 2001:2010)
#' # Interpolate the observations onto the model's grid. We use the method "nearest" and the getGrid function to ensure spatial consistency:
#' obs <- interpGridData(obs, new.grid = getGrid(prd), method = "nearest")
#' # Apply the bias correction method:
#' simBC <- isimip (obs, prd, sim, pr.threshold = 1) # ISI-MIP Method
#' par(mfrow = c(1,2))
#' plotMeanField(sim)
#' plotMeanField(simBC)
#' par(mfrow = c(1,1))
#' }

isimip <- function (obs, pred, sim, pr.threshold = 1) {
      datesObs <- as.POSIXct(obs$Dates$start, tz="GMT", format="%Y-%m-%d %H:%M:%S")
      datesObs<-cut(datesObs, "month")
      yearList<-unlist(strsplit(as.character(datesObs), "[-]"))
      months<-unique(unlist(strsplit(as.character(datesObs), "[-]"))[seq(2,length(yearList),3)])
      dimObs<-dim(obs$Data)
      dimMonthlyObs<-dim(obs$Data)
      dimMonthlyObs[which(grepl("^time",attr(obs$Data, "dimensions")))]<-length(unique(datesObs))
      dimPred<-dim(pred$Data)
      dimMonthlyPred<-dim(pred$Data)
      dimMonthlyPred[which(grepl("^time",attr(pred$Data, "dimensions")))]<-length(unique(datesObs))
      monthlyObs<- array(data = NA, dim = dimMonthlyObs)
      monthlyPred<- array(data = NA, dim = dimMonthlyPred)
      month2dayObs<- array(data = NA, dim = c(length(datesObs),1))
      obs.time.index <- grep("^time$", attr(obs$Data, "dimensions"))
      pred.time.index <- grep("^time$", attr(pred$Data, "dimensions"))
      for (i in 1:dimMonthlyObs[obs.time.index]){
            indMonth<-which(datesObs == unique(datesObs)[i])
            month2dayObs[indMonth]<-i
            indTimeObs <- rep(list(bquote()), length(dimObs))
            for (d in 1:length(dimObs)){
                  indTimeObs[[d]] <- 1:dimObs[d]
            }
            indTimeObs[[obs.time.index]] <- i
            indTimeObs<-as.matrix(expand.grid(indTimeObs))
            indTimeObs1 <- rep(list(bquote()), length(dimObs))
            indTimeObs1[[obs.time.index]] <- indMonth
            callObs <- as.call(c(list(as.name("["),quote(obs$Data)), indTimeObs1))
            monthlyObs[indTimeObs] <- apply(eval(callObs), FUN = mean, MARGIN = setdiff(1:length(dimObs),obs.time.index), na.rm = TRUE)
            indTimeObs <- rep(list(bquote()), length(dimPred))
            for (d in 1:length(dimPred)){
                  indTimeObs[[d]] <- 1:dimPred[d]
            }
            indTimeObs[[pred.time.index]] <- i
            indTimeObs<-as.matrix(expand.grid(indTimeObs))
            indTimeObs1 <- rep(list(bquote()), length(dimPred))
            indTimeObs1[[pred.time.index]] <- indMonth
            callObs <- as.call(c(list(as.name("["),quote(pred$Data)), indTimeObs1))
            monthlyPred[indTimeObs] <- apply(eval(callObs), FUN = mean, MARGIN = setdiff(1:length(dimPred),pred.time.index), na.rm = TRUE)
      }
      datesFor <- as.POSIXct(sim$Dates$start, tz="GMT", format="%Y-%m-%d %H:%M:%S")
      datesFor<-cut(datesFor, "month")
      yearListFor<-unlist(strsplit(as.character(datesFor), "[-]"))
      dimFor<-dim(sim$Data)
      dimMonthlyFor<-dim(sim$Data)
      dimMonthlyFor[which(grepl("^time",attr(sim$Data, "dimensions")))]<-length(unique(datesFor))
      monthlyFor<- array(data = NA, dim = dimMonthlyFor)
      month2dayFor<- array(data = NA, dim = c(length(datesFor),1))
      sim.time.index <- grep("^time$", attr(sim$Data, "dimensions"))
      for (i in 1:dimMonthlyFor[sim.time.index]){
            indMonth<-which(datesFor == unique(datesFor)[i])
            month2dayFor[indMonth]<-i
            indTimeObs <- rep(list(bquote()), length(dimFor))
            for (d in 1:length(dimFor)){
                  indTimeObs[[d]] <- 1:dimFor[d]
            }
            indTimeObs[[sim.time.index]] <- i
            indTimeObs<-as.matrix(expand.grid(indTimeObs))
            indTimeObs1 <- rep(list(bquote()), length(dimFor))
            indTimeObs1[[sim.time.index]] <- indMonth
            callObs <- as.call(c(list(as.name("["),quote(sim$Data)), indTimeObs1))
            monthlyFor[indTimeObs] <- apply(eval(callObs), FUN = mean, MARGIN = setdiff(1:length(dimFor),sim.time.index), na.rm = TRUE)
      }
      if (any(grepl(obs$Variable$varName,c("tas","mean temperature","tmean")))){
            # First Step: Monthly Correction
            dimAux<-dimPred
            dimAux[pred.time.index]<-length(months)
            monthlyCorrection<-array(data = NA, dim = dimAux)
            for (i in 1:length(months)){
                  indMonth <- which(grepl(paste("-",months[i],"-", sep=""),unique(datesObs)))
                  indTimeObs <- rep(list(bquote()), length(dimObs))
                  indTimeObs[[obs.time.index]] <- indMonth
                  callObs <- as.call(c(list(as.name("["),quote(monthlyObs)), indTimeObs))
                  auxObs <- apply(eval(callObs), FUN = mean, MARGIN = setdiff(1:length(dimObs),obs.time.index), na.rm = TRUE)
                  indTimeObs <- rep(list(bquote()), length(dimPred))
                  indTimeObs[[pred.time.index]] <- indMonth
                  callObs <- as.call(c(list(as.name("["),quote(monthlyPred)), indTimeObs))
                  auxPrd <- apply(eval(callObs), FUN = mean, MARGIN = setdiff(1:length(dimPred),pred.time.index), na.rm = TRUE)
                  indTimeObs <- rep(list(bquote()), length(dimPred))
                  for (d in 1:length(dimPred)){
                        indTimeObs[[d]] <- 1:dimPred[d]
                  }
                  indTimeObs[[pred.time.index]] <- i
                  indTimeObs<-as.matrix(expand.grid(indTimeObs))
                  auxObs <- array(auxObs, dim = c(dim(auxObs),dimPred[setdiff(1:length(dimPred),match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions")))]))
                  monthlyCorrection[indTimeObs]<-aperm(auxObs,match(attr(pred$Data,"dimensions")[setdiff(1:length(dimPred),pred.time.index)],c(attr(obs$Data,"dimensions")[setdiff(1:length(dimObs),obs.time.index)],attr(pred$Data,"dimensions")[setdiff(1:length(dimPred),match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions")))])))-auxPrd
            }
            indTimeObs <- rep(list(bquote()), length(dimPred))
            for (d in 1:length(dimPred)){
                  indTimeObs[[d]] <- 1:dimPred[d]
            }
            indTimeObs[[pred.time.index]] <- month2dayObs
            indTimeObs<-as.matrix(expand.grid(indTimeObs))
            pred$Data <- pred$Data-array(monthlyPred[indTimeObs], dim = dimPred)
            indTimeObs <- rep(list(bquote()), length(dimObs))
            for (d in 1:length(dimObs)){
                  indTimeObs[[d]] <- 1:dimObs[d]
            }
            indTimeObs[[obs.time.index]] <- month2dayObs
            indTimeObs<-as.matrix(expand.grid(indTimeObs))
            obs$Data <- obs$Data-array(monthlyObs[indTimeObs], dim = dimObs)
            indTimeObs <- rep(list(bquote()), length(dimFor))
            for (d in 1:length(dimFor)){
                  indTimeObs[[d]] <- 1:dimFor[d]
            }
            indTimeObs[[sim.time.index]] <- month2dayFor
            indTimeObs<-as.matrix(expand.grid(indTimeObs))
            sim$Data <- sim$Data-array(monthlyFor[indTimeObs], dim = dimFor)      
            indTimeObs <- rep(list(bquote()), length(dimObs))
            indTimePrd <- rep(list(bquote()), length(dimPred))
            indTimeSim <- rep(list(bquote()), length(dimFor))
            for (d in 1:length(dimObs)){
                  if (d!=obs.time.index){
                        indTimeObs[[d]] <- 1:dimObs[d]
                        indTimePrd[[match(attr(obs$Data, "dimensions")[d], attr(pred$Data, "dimensions"))]] <- 1:dimPred[match(attr(obs$Data, "dimensions")[d], attr(pred$Data, "dimensions"))]
                        indTimeSim[[match(attr(obs$Data, "dimensions")[d], attr(sim$Data, "dimensions"))]] <- 1:dimFor[match(attr(obs$Data, "dimensions")[d], attr(sim$Data, "dimensions"))]
                  }
            }
            dimDiff <- setdiff(1:length(dimPred),match(attr(obs$Data, "dimensions"), attr(pred$Data, "dimensions")))
            if (length(dimDiff)>0){
                  for (d in 1:length(dimDiff)){
                        indTimePrd[[dimDiff[d]]] <- 1:dimPred[dimDiff[d]]
                        indTimeSim[[dimDiff[d]]] <- 1:dimFor[dimDiff[d]]
                  }
            }
            indTimeObs[[obs.time.index]]<-NA
            indTimePrd[[pred.time.index]]<-NA
            indTimeSim[[sim.time.index]]<-NA
            indTimeObs<-as.matrix(expand.grid(indTimeObs))
            indTimePrd<-as.matrix(expand.grid(indTimePrd))
            indTimeSim<-as.matrix(expand.grid(indTimeSim))
            for (m in 1:length(months)){
                  indMonth<-which(grepl(paste("-",months[m],"-", sep=""),datesObs))
                  indMonthFor<-which(grepl(paste("-",months[m],"-", sep=""),datesFor))
                  for (j in 1:dim(indTimeObs)[1]){
                        indObs <- as.list(indTimeObs[j,])
                        indObs[[obs.time.index]]<-indMonth
                        indObs<-as.matrix(expand.grid(indObs))
                        if (any(!is.na(obs$Data[indObs]))){
                              indPrd <- rep(list(bquote()), length(dimPred))
                              indSim <- rep(list(bquote()), length(dimFor))
                              indPrd[match(attr(obs$Data, "dimensions"), attr(pred$Data, "dimensions"))] <- as.list(indTimeObs[j,])
                              indSim[match(attr(obs$Data, "dimensions"), attr(sim$Data, "dimensions"))] <- as.list(indTimeObs[j,])
                              if (length(dimDiff)>0){
                                    for (d in 1:length(dimDiff)){
                                          indPrd[[dimDiff[d]]] <- 1:dimPred[dimDiff[d]]
                                          indSim[[dimDiff[d]]] <- 1:dimFor[dimDiff[d]]
                                    }
                              }
                              indPrd<-as.matrix(expand.grid(indPrd))
                              indSim<-as.matrix(expand.grid(indSim))
                              for (i in 1:dim(indPrd)[1]){
                                    indPrd1 <- as.list(indPrd[i,])
                                    indPrd1[[pred.time.index]]<-indMonth
                                    indPrd1<-as.matrix(expand.grid(indPrd1))
                                    if (any(!is.na(pred$Data[indPrd1]))){
                                          indSim1 <- as.list(indSim[i,])
                                          indSim1[[sim.time.index]]<-indMonthFor
                                          indSim1<-as.matrix(expand.grid(indSim1))
                                          auxlmPrd <- pred$Data[indPrd1]
                                          auxlmObs <- obs$Data[indObs]
                                          indLM <- which(!is.na(auxlmPrd) & !is.na(auxlmObs))
                                          if (length(indLM)>0){
                                                lmAdjust<-lm(sort(auxlmPrd[indLM], decreasing = FALSE, na.last = NA) ~ sort(auxlmObs[indLM], decreasing = FALSE, na.last = NA)-1)
                                                sim$Data[indSim1]<-coef(lmAdjust)[1]*sim$Data[indSim1]
                                          }
                                    }
                              }
                        }
                  }
            }
            indTimeFor <- rep(list(bquote()), length(dimFor))
            for (d in 1:length(dimFor)){
                  indTimeFor[[d]] <- 1:dimFor[d]
            }
            indTimeFor[[sim.time.index]] <- month2dayFor
            indTimeFor<-as.matrix(expand.grid(indTimeFor))
            sim$Data <- sim$Data+array(monthlyFor[indTimeFor], dim = dimFor)      
            for (i in 1:length(months)){
                  indMonthFor<-which(grepl(paste("-",months[i],"-", sep=""),datesFor))
                  if (!is.null(indMonthFor)){
                        indCorrection <- rep(list(bquote()), length(dim(monthlyCorrection)))
                        for (d in 1:length(dimFor)){
                              indCorrection[[d]] <- 1:dim(monthlyCorrection)[d]
                        }
                        indCorrection[[pred.time.index]] <- i
                        indCorrection<-as.matrix(expand.grid(indCorrection))
                        for (j in 1:length(indMonthFor)){
                              indTimeFor <- rep(list(bquote()), length(dimFor))
                              for (d in 1:length(dimFor)){
                                    indTimeFor[[d]] <- 1:dimFor[d]
                              }
                              indTimeFor[[sim.time.index]] <- indMonthFor[j]
                              indTimeFor<-as.matrix(expand.grid(indTimeFor))
                              sim$Data[indTimeFor]<-sim$Data[indTimeFor]+monthlyCorrection[indCorrection]
                        }
                  }
            }
      }
      if (any(grepl(obs$Variable$varName,c("pr","tp","precipitation","precip")))){
            threshold<-pr.threshold
            if (length(threshold)==1){
                  threshold<-array(data = threshold, dim = 3)
            }
            # First Step: Monthly Correction
            dimAux<-dimPred
            dimAux[pred.time.index]<-length(months)
            monthlyCorrection<-array(data = NA, dim = dimAux)
            for (i in 1:length(months)){
                  indMonth<-which(grepl(paste("-",months[i],"-", sep=""),unique(datesObs)))
                  indTimeObs <- rep(list(bquote()), length(dimObs))
                  indTimeObs[[obs.time.index]] <- indMonth
                  callObs <- as.call(c(list(as.name("["),quote(monthlyObs)), indTimeObs))
                  auxObs <- apply(eval(callObs), FUN = sum, MARGIN = setdiff(1:length(dimObs),obs.time.index), na.rm = TRUE)
                  indTimeObs <- rep(list(bquote()), length(dimPred))
                  indTimeObs[[pred.time.index]] <- indMonth
                  callObs <- as.call(c(list(as.name("["),quote(monthlyPred)), indTimeObs))
                  auxPrd <- apply(eval(callObs), FUN = sum, MARGIN = setdiff(1:length(dimPred),pred.time.index), na.rm = TRUE)
                  indTimeObs <- rep(list(bquote()), length(dimPred))
                  for (d in 1:length(dimPred)){
                        indTimeObs[[d]] <- 1:dimPred[d]
                  }
                  indTimeObs[[pred.time.index]] <- i
                  indTimeObs<-as.matrix(expand.grid(indTimeObs))
                  auxObs <- array(auxObs, dim = c(dim(auxObs),dimPred[setdiff(1:length(dimPred),match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions")))]))
                  monthlyCorrection[indTimeObs]<-aperm(auxObs,match(dim(auxPrd),dim(auxObs)))/auxPrd
            }
            # Wet/dry days frequency correction:
            obs2perm <- match(attr(pred$Data, "dimensions"),attr(obs$Data, "dimensions"))[which(!is.na(match(attr(pred$Data, "dimensions"),attr(obs$Data, "dimensions"))))]
            auxPermObs <- 1:length(dimPred)
            auxPermObs[setdiff(1:length(dimPred),match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions")))] <- length(dimObs)+1:length(setdiff(1:length(dimPred),match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions"))))
            auxPermObs[match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions"))] <- obs2perm
            auxObs <- aperm(array(monthlyObs, dim = c(dim(monthlyObs),dimPred[setdiff(1:length(dimPred),match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions")))])),auxPermObs)
            nP <- apply(as.array(auxObs<=threshold[1] & !is.na(auxObs))*as.array(monthlyPred<=0 & !is.na(monthlyPred)), FUN = sum, MARGIN = setdiff(1:length(dimPred),pred.time.index))# Interpretamos la formula como interseccion
            Ndry <- array(data = NA, dim = dimObs[setdiff(1:length(dimObs), obs.time.index)])
            epsM <- array(data = NA, dim = dimPred[setdiff(1:length(dimPred), pred.time.index)])
            indTimePrd <- rep(list(bquote()), length(setdiff(1:length(dimPred),pred.time.index)))
            for (d in 1:length(setdiff(1:length(dimPred),pred.time.index))){
                  indTimePrd[[d]] <- 1:dimPred[setdiff(1:length(dimPred),pred.time.index)[d]]
            }
            indTimePrd<-as.matrix(expand.grid(indTimePrd))
            for (i in 1:dim(indTimePrd)[1]){
                  indTimePrdi <- rep(list(bquote()), length(dimPred))
                  indTimePrdi[[pred.time.index]] <- 1:dimMonthlyPred[pred.time.index]
                  indTimePrdi[setdiff(1:length(dimPred),pred.time.index)] <- as.list(indTimePrd[i,])
                  indTimePrdi<-as.matrix(expand.grid(indTimePrdi))
                  auxP <- sort(monthlyPred[indTimePrdi], decreasing = FALSE, na.last = NA)
                  indTimePrdi <- as.list(indTimePrd[i,])
                  indTimePrdi<-as.matrix(expand.grid(indTimePrdi))
                  epsM[indTimePrdi] <- auxP[min(c(nP[indTimePrdi]+1,dimPred[pred.time.index]), na.rm = TRUE)]
            }
            indTimePrd <- rep(list(bquote()), length(setdiff(1:length(dimObs),obs.time.index)))
            for (d in 1:length(setdiff(1:length(dimObs),obs.time.index))){
                  indTimePrd[[d]] <- 1:dimObs[setdiff(1:length(dimObs),obs.time.index)[d]]
            }
            indTimePrd<-as.matrix(expand.grid(indTimePrd))
            for (i in 1:dim(indTimePrd)[1]){
                  indTimePrdi <- rep(list(bquote()), length(dimObs))
                  indTimePrdi[[obs.time.index]] <- 1:dimMonthlyObs[obs.time.index]
                  indTimePrdi[setdiff(1:length(dimObs),obs.time.index)] <- as.list(indTimePrd[i,])
                  indTimePrdi<-as.matrix(expand.grid(indTimePrdi))
                  indMonth <- which(monthlyObs[indTimePrdi]>threshold[1]);
                  if (length(indMonth)>0){
                        indWetMonth <- which(month2dayObs%in%indMonth)
                        indTimePrdi <- rep(list(bquote()), length(dimObs))
                        indTimePrdi[[obs.time.index]] <- indWetMonth
                        indTimePrdi[setdiff(1:length(dimObs),obs.time.index)] <- as.list(indTimePrd[i,])
                        indTimePrdi<-as.matrix(expand.grid(indTimePrdi))
                        auxP <- obs$Data[indTimePrdi]
                        indTimePrdi <- as.list(indTimePrd[i,])
                        indTimePrdi<-as.matrix(expand.grid(indTimePrdi))
                        Ndry[indTimePrdi] <- length(which(auxP<threshold[3]))
                  }
            }
            epsD <- array(data = NA, dim = dimPred[setdiff(1:length(dimPred), pred.time.index)])
            indTimePrd <- rep(list(bquote()), length(setdiff(1:length(dimPred),pred.time.index)))
            for (d in 1:length(setdiff(1:length(dimPred),pred.time.index))){
                  indTimePrd[[d]] <- 1:dimPred[setdiff(1:length(dimPred),pred.time.index)[d]]
            }
            indTimePrd<-as.matrix(expand.grid(indTimePrd))
            for (i in 1:dim(indTimePrd)[1]){
                  indTimePrdi <- rep(list(bquote()), length(dimPred))
                  indTimePrdi[[pred.time.index]] <- 1:dimMonthlyPred[pred.time.index]
                  indTimePrdi[setdiff(1:length(dimPred),pred.time.index)] <- as.list(indTimePrd[i,])
                  indTimePrdi<-as.matrix(expand.grid(indTimePrdi))
                  auxP <- monthlyPred[indTimePrdi]
                  indTimePrdi <- as.list(indTimePrd[i,])
                  indTimePrdi<-as.matrix(expand.grid(indTimePrdi))
                  indMonth <- which(auxP>epsM[indTimePrdi])
                  if (length(indMonth)>0){
                        indWetMonth <- which(month2dayObs%in%indMonth)
                        indWet <- rep(list(bquote()), length(dimPred))
                        indWet[[pred.time.index]] <- indWetMonth
                        indWet[setdiff(1:length(dimPred),pred.time.index)] <- as.list(indTimePrd[i,])
                        indWet <- as.matrix(expand.grid(indWet))
                        indNdry <- as.list(indTimePrd[i,which(match(attr(pred$Data, "dimensions")[setdiff(1:length(dimPred),pred.time.index)],attr(obs$Data, "dimensions")[setdiff(1:length(dimObs),obs.time.index)], nomatch = 0)>0)])
                        indNdry <- as.matrix(expand.grid(indNdry))
                        auxP <- sort(pred$Data[indWet], decreasing = FALSE, na.last = NA)
                        if (is.na(Ndry[indNdry]) | Ndry[indNdry]==length(indWetMonth)){
                              Ndry[indNdry] <- length(indWetMonth)
                              epsD[indTimePrdi] <- max(pred$Data[indWet[which(pred$Data[indWet]<=auxP[Ndry[indNdry]]),]], na.rm = TRUE)
                        }else{
                              epsD[indTimePrdi] <- 0.5*max(pred$Data[indWet[which(pred$Data[indWet]<=auxP[Ndry[indNdry]]),]], na.rm = TRUE)+0.5*min(pred$Data[indWet[which(pred$Data[indWet]>auxP[Ndry[indNdry]]),]], na.rm = TRUE)
                        }
                        for (j in 1:length(indMonth)){
                              indWetMonth <- which(month2dayObs==indMonth[j])
                              indWet <- rep(list(bquote()), length(dimPred))
                              indWet[[pred.time.index]] <- indWetMonth
                              indWet[setdiff(1:length(dimPred),pred.time.index)] <- as.list(indTimePrd[i,])
                              indWet <- as.matrix(expand.grid(indWet))
                              indDryMonth <- indWetMonth[which(pred$Data[indWet]<=epsD[indTimePrdi])]
                              if (length(indDryMonth)>0){
                                    indDry <- rep(list(bquote()), length(dimPred))
                                    indDry[[pred.time.index]] <- indDryMonth
                                    indDry[setdiff(1:length(dimPred),pred.time.index)] <- as.list(indTimePrd[i,])
                                    indDry <- as.matrix(expand.grid(indDry))
                                    mi <- sum(pred$Data[indDry], na.rm= TRUE)/(length(indWetMonth)-length(indDryMonth))
                                    pred$Data[indWet[which(pred$Data[indWet]>epsD[indTimePrdi]),]] <- pred$Data[indWet[which(pred$Data[indWet]>epsD[indTimePrdi]),]]+mi
                                    pred$Data[indDry] <- 0
                              }
                        }
                  }
            }
            indTimePrd <- rep(list(bquote()), length(setdiff(1:length(dimFor),sim.time.index)))
            for (d in 1:length(setdiff(1:length(dimFor),sim.time.index))){
                  indTimePrd[[d]] <- 1:dimFor[setdiff(1:length(dimFor),sim.time.index)[d]]
            }
            indTimePrd<-as.matrix(expand.grid(indTimePrd))
            for (i in 1:dim(indTimePrd)[1]){
                  indTimePrdi <- rep(list(bquote()), length(dimFor))
                  indTimePrdi[[sim.time.index]] <- 1:dimMonthlyPred[sim.time.index]
                  indTimePrdi[setdiff(1:length(dimFor),sim.time.index)] <- as.list(indTimePrd[i,])
                  indTimePrdi<-as.matrix(expand.grid(indTimePrdi))
                  auxP <- monthlyFor[indTimePrdi]
                  indTimePrdi <- as.list(indTimePrd[i,])
                  indTimePrdi<-as.matrix(expand.grid(indTimePrdi))
                  indMonth <- which(auxP>epsM[indTimePrdi])
                  if (length(indMonth)>0){
                        for (j in 1:length(indMonth)){
                              indWetMonth <- which(month2dayFor==indMonth[j])
                              indWet <- rep(list(bquote()), length(dimFor))
                              indWet[[sim.time.index]] <- indWetMonth
                              indWet[setdiff(1:length(dimFor),sim.time.index)] <- as.list(indTimePrd[i,])
                              indWet <- as.matrix(expand.grid(indWet))
                              indDryMonth <- indWetMonth[which(sim$Data[indWet]<=epsD[indTimePrdi])]
                              if (length(indDryMonth)>0){
                                    indDry <- rep(list(bquote()), length(dimFor))
                                    indDry[[sim.time.index]] <- indDryMonth
                                    indDry[setdiff(1:length(dimFor),sim.time.index)] <- as.list(indTimePrd[i,])
                                    indDry <- as.matrix(expand.grid(indDry))
                                    mi <- sum(sim$Data[indDry], na.rm= TRUE)/(length(indWetMonth)-length(indDryMonth))
                                    sim$Data[indWet[which(sim$Data[indWet]>epsD[indTimePrdi]),]] <- sim$Data[indWet[which(sim$Data[indWet]>epsD[indTimePrdi]),]]+mi
                                    sim$Data[indDry] <- 0
                              }
                        }
                  }
            }
            indTimeObs <- rep(list(bquote()), length(dimObs))
            for (d in 1:length(dimObs)){
                  indTimeObs[[d]] <- 1:dimObs[d]
            }
            indTimeObs[[obs.time.index]] <- month2dayObs
            indTimeObs<-as.matrix(expand.grid(indTimeObs))
            indTimePrd <- rep(list(bquote()), length(dimPred))
            for (d in 1:length(dimPred)){
                  indTimePrd[[d]] <- 1:dimPred[d]
            }
            indTimePrd[[pred.time.index]] <- month2dayObs
            indTimePrd<-as.matrix(expand.grid(indTimePrd))
            indTimeFor <- rep(list(bquote()), length(dimFor))
            for (d in 1:length(dimFor)){
                  indTimeFor[[d]] <- 1:dimFor[d]
            }
            indTimeFor[[sim.time.index]] <- month2dayFor
            indTimeFor<-as.matrix(expand.grid(indTimeFor))
            auxPerm <- 1:length(dimPred)
            auxPerm[pred.time.index] <- length(dim(epsD))+1
            auxPerm[setdiff(1:length(dimPred),pred.time.index)] <- 1:length(dim(epsD))
            auxPermObs <- 1:length(dimPred)
            auxPermObs[setdiff(1:length(dimPred),match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions")))] <- length(dimObs)+1:length(setdiff(1:length(dimPred),match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions"))))
            auxPermObs[match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions"))] <- obs2perm
            wetDaysObs <- as.array(!is.na(pred$Data) & pred$Data>aperm(array(epsD, dim = c(dim(epsD),dimPred[pred.time.index])), auxPerm))*as.array(monthlyPred[indTimePrd]>aperm(array(epsM, dim = c(dim(epsD),dimPred[pred.time.index])), auxPerm))*aperm(array(as.array(obs$Data> threshold[3] & !is.na(obs$Data)), dim = c(dimPred[match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions"))],dimPred[setdiff(1:length(dimPred),match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions")))])),auxPermObs)*aperm(array(as.array(array(monthlyObs[indTimeObs], dim = dimObs)> threshold[1] & !is.na(as.array(array(monthlyObs[indTimeObs], dim = dimObs)))), dim = c(dimPred[match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions"))],dimPred[setdiff(1:length(dimPred),match(attr(obs$Data,"dimensions"),attr(pred$Data,"dimensions")))])),auxPermObs)
            wetDaysFor <- as.array(!is.na(sim$Data) & sim$Data>aperm(array(epsD, dim = c(dim(epsD),dimFor[sim.time.index])), auxPerm))*as.array(monthlyFor[indTimeFor]>aperm(array(epsM, dim = c(dim(epsD),dimFor[sim.time.index])), auxPerm))
            indTime <- rep(list(bquote()), length(setdiff(1:length(dimFor),sim.time.index)))
            for (d in 1:length(setdiff(1:length(dimFor),sim.time.index))){
                  indTime[[d]] <- 1:dimFor[setdiff(1:length(dimFor),sim.time.index)[d]]
            }
            indTime<-as.matrix(expand.grid(indTime))
            for (i in 1:dim(indTime)[1]){
                  indTimeFor <- rep(list(bquote()), length(dimFor))
                  indTimeFor[[sim.time.index]] <- month2dayFor
                  indTimeFor[setdiff(1:length(dimFor),sim.time.index)] <- as.list(indTime[i,])
                  indTimeFor<-as.matrix(expand.grid(indTimeFor))
                  auxMonthly <- monthlyFor[indTimeFor]
                  indAuxMonthly <- which(auxMonthly > epsM[as.matrix(expand.grid(as.list(indTime[i,])))])
                  indTimeFor <- rep(list(bquote()), length(dimFor))
                  indTimeFor[[sim.time.index]] <- indAuxMonthly
                  indTimeFor[setdiff(1:length(dimFor),sim.time.index)] <- as.list(indTime[i,])
                  indTimeFor<-as.matrix(expand.grid(indTimeFor))
                  sim$Data[indTimeFor] <- sim$Data[indTimeFor]/auxMonthly[indAuxMonthly]
                  indTimeFor <- rep(list(bquote()), length(dimPred))
                  indTimeFor[[pred.time.index]] <- 1:dimPred[pred.time.index]
                  indTimeFor[setdiff(1:length(dimPred),pred.time.index)] <- as.list(indTime[i,])
                  indTimeFor<-as.matrix(expand.grid(indTimeFor))
                  indAuxMonthly <- which(wetDaysObs[indTimeFor]==1)
                  indTimeFor <- rep(list(bquote()), length(dimPred))
                  indTimeFor[[pred.time.index]] <- indAuxMonthly
                  indTimeFor[setdiff(1:length(dimPred),pred.time.index)] <- as.list(indTime[i,])
                  indTimeFor<-as.matrix(expand.grid(indTimeFor))
                  indTimeObs <- rep(list(bquote()), length(dimPred))
                  indTimeObs[[pred.time.index]] <- month2dayObs[indAuxMonthly]
                  indTimeObs[setdiff(1:length(dimPred),pred.time.index)] <- as.list(indTime[i,])
                  indTimeObs<-as.matrix(expand.grid(indTimeObs))
                  pred$Data[indTimeFor] <- pred$Data[indTimeFor]/monthlyPred[indTimeObs]
            }
            indTime <- rep(list(bquote()), length(setdiff(1:length(dimObs),obs.time.index)))
            for (d in 1:length(setdiff(1:length(dimObs),obs.time.index))){
                  indTime[[d]] <- 1:dimObs[setdiff(1:length(dimObs),obs.time.index)[d]]
            }
            indTime<-as.matrix(expand.grid(indTime))
            indTimePrd <- rep(list(bquote()), length(setdiff(1:length(dimPred),pred.time.index)))
            for (d in 1:length(setdiff(1:length(dimPred),pred.time.index))){
                  if (is.na(match(attr(pred$Data, "dimensions")[setdiff(1:length(dimPred),pred.time.index)[d]],attr(obs$Data, "dimensions")))){
                        indTimePrd[[d]] <- 1
                  }else{
                        indTimePrd[[d]] <- 1:dimPred[setdiff(1:length(dimPred),pred.time.index)[d]]
                  }
            }
            indTimePrd<-as.matrix(expand.grid(indTimePrd))
            for (i in 1:dim(indTime)[1]){
                  indTimeFor <- rep(list(bquote()), length(dimPred))
                  indTimeFor[[pred.time.index]] <- 1:dimPred[pred.time.index]
                  indTimeFor[setdiff(1:length(dimPred),pred.time.index)] <- as.list(indTimePrd[i,])
                  indTimeFor<-as.matrix(expand.grid(indTimeFor))
                  indAuxMonthly <- which(wetDaysObs[indTimeFor]==1)
                  indTimeFor <- rep(list(bquote()), length(dimObs))
                  indTimeFor[[obs.time.index]] <- indAuxMonthly
                  indTimeFor[setdiff(1:length(dimObs),obs.time.index)] <- as.list(indTime[i,])
                  indTimeFor<-as.matrix(expand.grid(indTimeFor))
                  indTimeObs <- rep(list(bquote()), length(dimObs))
                  indTimeObs[[obs.time.index]] <- month2dayObs[indAuxMonthly]
                  indTimeObs[setdiff(1:length(dimObs),obs.time.index)] <- as.list(indTime[i,])
                  indTimeObs<-as.matrix(expand.grid(indTimeObs))
                  obs$Data[indTimeFor] <- obs$Data[indTimeFor]/monthlyObs[indTimeObs]
            }
            indEstObs <- rep(list(bquote()), length(setdiff(1:length(dimObs),obs.time.index)))
            for (d in 1:length(setdiff(1:length(dimObs),obs.time.index))){
                  indEstObs[[d]] <- 1:dimObs[setdiff(1:length(dimObs),obs.time.index)[d]]
            }
            indEstObs<-as.matrix(expand.grid(indEstObs))
            indEstPrd <- rep(list(bquote()), length(setdiff(1:length(dimPred),pred.time.index)))
            for (d in 1:length(setdiff(1:length(dimPred),pred.time.index))){
                  indEstPrd[[d]] <- 1:dimPred[setdiff(1:length(dimPred),pred.time.index)[d]]
            }
            indEstPrd<-as.matrix(expand.grid(indEstPrd))
            for (i in 1:dim(indEstPrd)[1]){
                  indTime <- rep(list(bquote()), length(dimPred))
                  indTime[[pred.time.index]] <- 1:dimPred[pred.time.index]
                  indTime[setdiff(1:length(dimPred),pred.time.index)] <- as.list(indEstPrd[i,])
                  indTime<-as.matrix(expand.grid(indTime))
                  indTimeSim <- rep(list(bquote()), length(dimFor))
                  indTimeSim[[sim.time.index]] <- 1:dimFor[sim.time.index]
                  indTimeSim[setdiff(1:length(dimFor),sim.time.index)] <- as.list(indEstPrd[i,])
                  indTimeSim<-as.matrix(expand.grid(indTimeSim))
                  for (j in 1:length(months)){
                        indMonth <- which(unlist(strsplit(as.character(datesObs), "[-]"))[seq(2,length(yearList),3)] == months[j] & wetDaysObs[indTime] == 1)
                        indMonthFor <- which(unlist(strsplit(as.character(datesFor), "[-]"))[seq(2,length(yearListFor),3)] == months[j] & wetDaysFor[indTimeSim] == 1)
                        if (length(indMonth)>80 & length(indMonthFor)>0){
                              indTimePred <- rep(list(bquote()), length(dimPred))
                              indTimePred[[pred.time.index]] <- indMonth
                              indTimePred[setdiff(1:length(dimPred),pred.time.index)] <- as.list(indEstPrd[i,])
                              indTimePred<-as.matrix(expand.grid(indTimePred))
                              indTimeObs <- rep(list(bquote()), length(dimObs))
                              indTimeObs[[obs.time.index]] <- indMonth
                              indTimeObs[match(attr(pred$Data, "dimensions")[setdiff(1:length(dimPred),pred.time.index)], attr(obs$Data, "dimensions"))[!is.na(match(attr(pred$Data, "dimensions")[setdiff(1:length(dimPred),pred.time.index)], attr(obs$Data, "dimensions")))]] <- as.list(indEstPrd[i,match(attr(obs$Data, "dimensions")[setdiff(1:length(dimObs),obs.time.index)], attr(pred$Data, "dimensions")[setdiff(1:length(dimPred),pred.time.index)])])
                              indTimeObs<-as.matrix(expand.grid(indTimeObs))
                              indTimeFor <- rep(list(bquote()), length(dimFor))
                              indTimeFor[[sim.time.index]] <- indMonthFor
                              indTimeFor[setdiff(1:length(dimFor),sim.time.index)] <- as.list(indEstPrd[i,])
                              indTimeFor<-as.matrix(expand.grid(indTimeFor))
                              if (min(c(mean(pred$Data[indTimePred], na.rm = TRUE), mean(obs$Data[indTimeObs], na.rm = TRUE)), na.rm = TRUE)>0.01){
                                    auxO <- sort(obs$Data[indTimeObs], decreasing = FALSE, na.last = NA)
                                    auxP <- sort(pred$Data[indTimePred], decreasing = FALSE, na.last = NA)
                                    tryCatch({
                                          fit <- nls(auxO ~ (q1+q2*(auxP-min(auxP, na.rm = TRUE)))*(1-exp(-(auxP-min(auxP, na.rm = TRUE))/q3)),start=list(q1=runif(1, min = 0, max = 1),q2=runif(1, min = 0, max = 1),q3=runif(1, min = 0, max = 1)*(max(auxP, na.rm = TRUE)-min(auxP, na.rm = TRUE))), na.action = na.omit)
                                          if (sum(resid(fit)^2)<=1e-8 & length(auxP)>3){
                                                auxS <- abs(predict(fit,sim$Data[indTimeFor]))
                                          }else{
                                                lmAdjust<-lm(auxP ~ auxO)
                                                auxS <- abs(coef(lmAdjust)[1]+coef(lmAdjust)[2]*as.numeric(sim$Data[indTimeFor]))
                                          }
                                    }, error = function(e) e, finally = {
                                          lmAdjust<-lm(auxP ~ auxO)
                                          auxS <- abs(coef(lmAdjust)[1]+coef(lmAdjust)[2]*as.numeric(sim$Data[indTimeFor]))
                                    })
                                    sim$Data[indTimeFor] <- auxS
                              }
                        }
                  }
            }
            indTimeFor <- rep(list(bquote()), length(dimFor))
            for (d in 1:length(dimFor)){
                  indTimeFor[[d]] <- 1:dimFor[d]
            }
            indTimeFor[[sim.time.index]] <- month2dayFor
            indTimeFor<-as.matrix(expand.grid(indTimeFor))
            sim$Data <- sim$Data*array(monthlyFor[indTimeFor], dim = dimFor)      
            for (i in 1:length(months)){
                  indMonthFor<-which(grepl(paste("-",months[i],"-", sep=""),datesFor))
                  if (!is.null(indMonthFor)){
                        indCorrection <- rep(list(bquote()), length(dim(monthlyCorrection)))
                        for (d in 1:length(dimFor)){
                              indCorrection[[d]] <- 1:dim(monthlyCorrection)[d]
                        }
                        indCorrection[[pred.time.index]] <- i
                        indCorrection<-as.matrix(expand.grid(indCorrection))
                        for (j in 1:length(indMonthFor)){
                              indTimeFor <- rep(list(bquote()), length(dimFor))
                              for (d in 1:length(dimFor)){
                                    indTimeFor[[d]] <- 1:dimFor[d]
                              }
                              indTimeFor[[sim.time.index]] <- indMonthFor[j]
                              indTimeFor<-as.matrix(expand.grid(indTimeFor))
                              sim$Data[indTimeFor]<-sim$Data[indTimeFor]*monthlyCorrection[indCorrection]
                        }
                  }
            }
            attr(sim$Data, "threshold") <- threshold
      }
      attr(sim$Data, "correction") <- "ISI-MIP"
      return(sim)
}
