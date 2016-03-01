#' @title Bias correction ISI-MIP method
#' @description Implementation of the ISI-MIP methodology
#'
#' @template templateObsPredSim
#' @param threshold The minimum value that is considered as a non-occurrence (e.g. precipitation). Default to 1.
#' @param type Type of bias correction approach, either multiplicative (e.g. precipitation, \code{type = "multiplicative"})
#'  or additive (e.g. temperature, \code{type = "additive"}, the default).
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
#' @return A calibrated object of the same spatio-temporal extent of the input grid
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
#' obs <- loadStationData(dataset = value, var="precip", lonLim = c(-12,10), latLim = c(33,47),
#'                        season=c(12,1,2), years = 1991:2000)
#' prd <- loadGridData(ncep, var = "tp", lonLim = c(-12,10), latLim = c(33,47),
#'                     season = c(12,1,2), years = 1991:2000)
#' sim <- loadGridData(ncep, var = "tp", lonLim = c(-12,10), latLim = c(33,47),
#'                     season = c(12,1,2), years = 2001:2010)
#' # Interpolate the observations onto the model's grid. We use the method "nearest" 
#' # and the getGrid function to ensure spatial consistency:
#' obs <- interpData(obs, new.Coordinates = getGrid(prd), method = "nearest")
#' # Apply the bias correction method:
#' simBC <- isimip (obs, prd, sim, threshold = 1) # ISI-MIP Method
#' par(mfrow = c(1,2))
#' plotMeanGrid(sim)
#' plotMeanGrid(simBC)
#' par(mfrow = c(1,1))
#' }

isimip <- function (obs, pred, sim, threshold = 1, type = c("additive", "multiplicative")) {
      
      if (is.null(obs$Dates$start)){
            datesObs <- as.POSIXct(obs$Dates[[1]]$start, tz="GMT", format="%Y-%m-%d %H:%M:%S")
            multiField <- TRUE
            obs.var.index <- grep("^var$", attr(obs$Data, "dimensions"))
            pred.var.index <- grep("^var$", attr(pred$Data, "dimensions"))
            sim.var.index <- grep("^var$", attr(sim$Data, "dimensions"))
      }else{
            datesObs <- as.POSIXct(obs$Dates$start, tz="GMT", format="%Y-%m-%d %H:%M:%S")
            multiField <- FALSE
      }
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
      
      if (is.null(sim$Dates$start)){
            datesFor <- as.POSIXct(sim$Dates[[1]]$start, tz="GMT", format="%Y-%m-%d %H:%M:%S")            
      }else{
            datesFor <- as.POSIXct(sim$Dates$start, tz="GMT", format="%Y-%m-%d %H:%M:%S")            
      }
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
      if ((any(match(obs$Variable$varName,c("tas","mean temperature","tmean"))) | (type == "additive")) & !multiField){
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
      if (((any(match(obs$Variable$varName,c("radiation","pressure","wind","windspeed","humidity","specific humidity","radiacion","presion","viento","humedad","humedad especifica","rss","rsds","rls","rlds","ps","slp","wss","huss","hus")))) | (type == "multiplicative")) & !multiField){
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
            epsD <- array(data = 0, dim = dimPred[setdiff(1:length(dimPred), pred.time.index)])
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
      if (!is.na(any(match(obs$Variable$varName,c("maximum temperature","temperatura maxima","tasmax","tmax","minimum temperature","temperatura minima","tasmin","tmin"))))) {
            indTas <- which(!is.na(match(obs$Variable$varName,c("tas","temperatura media","mean temperature","tmean"))))
            if (length(indTas) == 0){
                  stop("Mean temperature is needed to correct the Maximum and Minimum Temperatures")
            }else{
                  if (!is.na(any(match(obs$Variable$varName,c("maximum temperature","temperatura maxima","tasmax","tmax"))))) {
                        indTasMax <- which(!is.na(match(obs$Variable$varName,c("maximum temperature","temperatura maxima","tasmax","tmax"))))
                  }
                  if (!is.na(any(match(obs$Variable$varName,c("minimum temperature","temperatura minima","tasmin","tmin"))))) {
                        indTasMax <- which(!is.na(match(obs$Variable$varName,c("minimum temperature","temperatura minima","tasmin","tmin"))))
                  }
                  obsTas <- obs
                  indTasAux <- rep(list(bquote()), length(dimObs))
                  indTasAux[[obs.var.index]] <- indTas
                  callTas <- as.call(c(list(as.name("["),quote(obs$Data)), indTasAux))
                  obsTas$Data <- eval(callTas)
                  attr(obsTas$Data, "dimensions") <- attributes(obs$Data)$dimensions[setdiff(1:length(attributes(obs$Data)$dimensions),obs.var.index)]
                  obsTas$Variable$varName <- obs$Variable$varName[indTas]
                  obsTas$Dates$start <- obs$Dates[[indTas]]$start
                  obsTas$Dates$end <- obs$Dates[[indTas]]$end
                  prdTas <- pred
                  indTasAux <- rep(list(bquote()), length(dimPred))
                  indTasAux[[pred.var.index]] <- indTas
                  callTas <- as.call(c(list(as.name("["),quote(pred$Data)), indTasAux))
                  prdTas$Data <- eval(callTas)
                  attr(prdTas$Data, "dimensions") <- attributes(pred$Data)$dimensions[setdiff(1:length(attributes(pred$Data)$dimensions),pred.var.index)]
                  prdTas$Variable$varName <- pred$Variable$varName[indTas]
                  prdTas$Dates$start <- pred$Dates[[indTas]]$start
                  prdTas$Dates$end <- pred$Dates[[indTas]]$end
                  simTas <- sim
                  indTasAux <- rep(list(bquote()), length(dimFor))
                  indTasAux[[sim.var.index]] <- indTas
                  callTas <- as.call(c(list(as.name("["),quote(sim$Data)), indTasAux))
                  simTas$Data <- eval(callTas)
                  attr(simTas$Data, "dimensions") <- attributes(sim$Data)$dimensions[setdiff(1:length(attributes(sim$Data)$dimensions),sim.var.index)]
                  simTas$Variable$varName <- sim$Variable$varName[indTas]
                  simTas$Dates$start <- sim$Dates[[indTas]]$start
                  simTas$Dates$end <- sim$Dates[[indTas]]$end
                  simTas <- isimip(obsTas, prdTas, simTas, type = "additive")# isimip
                  indTasFor <- rep(list(bquote()), length(dimFor))
                  for (d in 1:length(dimFor)){
                        indTasFor[[d]] <- 1:dimFor[d]
                  }
                  daysObs <- matrix(data = c(as.numeric(unlist(strsplit(as.character(datesObs), "[-]"))[seq(2,length(yearList),3)]),as.numeric(unlist(strsplit(as.character(datesObs), "[-]"))[seq(3,length(yearList),3)])), nrow = length(yearList)/3, ncol = 2)
                  daysObsU <- unique(daysObs, MARGIN = 1)
                  ndaysObs <- dim(daysObsU)[1]
                  daysFor <- matrix(data = c(as.numeric(unlist(strsplit(as.character(datesFor), "[-]"))[seq(2,length(yearListFor),3)]),as.numeric(unlist(strsplit(as.character(datesFor), "[-]"))[seq(3,length(yearListFor),3)])), nrow = length(yearListFor)/3, ncol = 2)
                  for (nd in 1:ndaysObs){
                        indDaysObs <- which(abs(daysObs[,1]-daysObsU[nd,1])+abs(daysObs[,2]-daysObsU[nd,2])==0)
                        if (!is.null(indDaysObs)){
                              indDaysFor <- which(abs(daysFor[,1]-daysObsU[nd,1])+abs(daysFor[,2]-daysObsU[nd,2])==0)
                              if (!is.null(indDaysFor)){
                                    indTasObsAux1 <- rep(list(bquote()), length(dimObs))
                                    indTasObsAux2 <- rep(list(bquote()), length(dimObs))
                                    indTasObsAux1[[obs.time.index]] <- indDaysObs
                                    indTasObsAux2[[obs.time.index]] <- indDaysObs
                                    indTasObsAux1[[obs.var.index]] <- indTas
                                    indTasObsAux2[[obs.var.index]] <- indTasMax
                                    callObs1 <- as.call(c(list(as.name("["),quote(obs$Data)), indTasObsAux1))
                                    callObs2 <- as.call(c(list(as.name("["),quote(obs$Data)), indTasObsAux2))
                                    indTasPrdAux1 <- rep(list(bquote()), length(dimPred))
                                    indTasPrdAux2 <- rep(list(bquote()), length(dimPred))
                                    indTasPrdAux1[[pred.time.index]] <- indDaysObs
                                    indTasPrdAux2[[pred.time.index]] <- indDaysObs
                                    indTasPrdAux1[[pred.var.index]] <- indTas
                                    indTasPrdAux2[[pred.var.index]] <- indTasMax
                                    callPrd1 <- as.call(c(list(as.name("["),quote(pred$Data)), indTasPrdAux1))
                                    callPrd2 <- as.call(c(list(as.name("["),quote(pred$Data)), indTasPrdAux2))
                                    MARGIN <- setdiff(1:length(attr(pred$Data, "dimensions")[setdiff(1:length(dimObs),obs.var.index)]),grep("^time$", attr(pred$Data, "dimensions")[setdiff(1:length(dimObs),obs.var.index)]))
                                    k <- apply(eval(callObs2)-eval(callObs1), FUN = sum, MARGIN = MARGIN, na.rm = TRUE)/
                                          apply(eval(callPrd2)-eval(callPrd1), FUN = sum, MARGIN = MARGIN, na.rm = TRUE)
                                    indTasFor[[sim.time.index]] <- indDaysFor
                                    indTasFor[[sim.var.index]] <- indTasMax
                                    indTasForAux <- as.matrix(expand.grid(indTasFor))
                                    indTasForAux1 <- rep(list(bquote()), length(dimFor))
                                    indTasForAux2 <- rep(list(bquote()), length(dimFor))
                                    indTasForAux3 <- rep(list(bquote()), length(dim(simTas$Data)))
                                    indTasForAux1[[sim.time.index]] <- indDaysFor
                                    indTasForAux2[[sim.time.index]] <- indDaysFor
                                    indTasForAux3[[grep("^time$", attr(simTas$Data, "dimensions"))]] <- indDaysFor
                                    indTasForAux1[[sim.var.index]] <- indTas
                                    indTasForAux2[[sim.var.index]] <- indTasMax
                                    callSim1 <- as.call(c(list(as.name("["),quote(sim$Data)), indTasForAux1))
                                    callSim2 <- as.call(c(list(as.name("["),quote(sim$Data)), indTasForAux2))
                                    callSim3 <- as.call(c(list(as.name("["),quote(simTas$Data)), indTasForAux3))
                                    k <- array(data = k, dim = c(dim(simTas$Data)[MARGIN],length(indDaysFor)))
                                    k <- aperm(k, perm = c(length(dim(k)), 1:(length(dim(k))-1)))
                                    sim$Data[indTasForAux] <- k*(eval(callSim2)-eval(callSim1))+eval(callSim3)
                              }
                        }
                  }
            }
      }
      #       case {'uas';'vas';'ua';'va';'eastward wind component';'northward wind component'},
      #       if isempty(Ws),error('Wind speed is necessary for the correction of the eastward and northward wind component');end
      #       wsC=isimip(Ws.O,Ws.P,Ws.F,'variable','windspeed','datesobs',datesObs,'datesfor',datesFor);
      #       indC=find(~isnan(Ws.F) & Ws.F>0);F(indC)=(F(indC).*wsC(indC))./Ws.F(indC);
      if (!is.na(any(match(obs$Variable$varName,c("uas","vas","ua","va","eastward wind component","northward wind component"))))) {
            indTas <- which(!is.na(match(obs$Variable$varName,c("wind","windspeed","viento","wss"))))
            if (length(indTas) == 0){
                  stop("Wind speed is needed to correct eastward and northward wind components")
            }else{
                  if (!is.na(any(match(obs$Variable$varName,c("uas","ua","eastward wind component"))))) {
                        indTasMax <- which(!is.na(match(obs$Variable$varName,c("uas","ua","eastward wind component"))))
                  }
                  if (!is.na(any(match(obs$Variable$varName,c("vas","va","northward wind component"))))) {
                        indTasMax <- which(!is.na(match(obs$Variable$varName,c("vas","va","northward wind component"))))
                  }
                  obsTas <- obs
                  indTasAux <- rep(list(bquote()), length(dimObs))
                  indTasAux[[obs.var.index]] <- indTas
                  callTas <- as.call(c(list(as.name("["),quote(obs$Data)), indTasAux))
                  obsTas$Data <- eval(callTas)
                  attr(obsTas$Data, "dimensions") <- attributes(obs$Data)$dimensions[setdiff(1:length(attributes(obs$Data)$dimensions),obs.var.index)]
                  obsTas$Variable$varName <- obs$Variable$varName[indTas]
                  obsTas$Dates$start <- obs$Dates[[indTas]]$start
                  obsTas$Dates$end <- obs$Dates[[indTas]]$end
                  prdTas <- pred
                  indTasAux <- rep(list(bquote()), length(dimPred))
                  indTasAux[[pred.var.index]] <- indTas
                  callTas <- as.call(c(list(as.name("["),quote(pred$Data)), indTasAux))
                  prdTas$Data <- eval(callTas)
                  attr(prdTas$Data, "dimensions") <- attributes(pred$Data)$dimensions[setdiff(1:length(attributes(pred$Data)$dimensions),pred.var.index)]
                  prdTas$Variable$varName <- pred$Variable$varName[indTas]
                  prdTas$Dates$start <- pred$Dates[[indTas]]$start
                  prdTas$Dates$end <- pred$Dates[[indTas]]$end
                  simTas <- sim
                  indTasAux <- rep(list(bquote()), length(dimFor))
                  indTasAux[[sim.var.index]] <- indTas
                  callTas <- as.call(c(list(as.name("["),quote(sim$Data)), indTasAux))
                  simTas$Data <- eval(callTas)
                  attr(simTas$Data, "dimensions") <- attributes(sim$Data)$dimensions[setdiff(1:length(attributes(sim$Data)$dimensions),sim.var.index)]
                  simTas$Variable$varName <- sim$Variable$varName[indTas]
                  simTas$Dates$start <- sim$Dates[[indTas]]$start
                  simTas$Dates$end <- sim$Dates[[indTas]]$end
                  simTas <- isimip(obsTas, prdTas, simTas, type = "multiplicative")# isimip
                  indTasFor <- rep(list(bquote()), length(dimFor))
                  for (d in 1:length(dimFor)){
                        indTasFor[[d]] <- 1:dimFor[d]
                  }
                  indTasForAux <- as.matrix(expand.grid(indTasFor))
                  indTasForAux1 <- rep(list(bquote()), length(dimFor))
                  indTasForAux2 <- rep(list(bquote()), length(dimFor))
                  indTasForAux1[[sim.var.index]] <- indTas
                  indTasForAux2[[sim.var.index]] <- indTasMax
                  callSim1 <- as.call(c(list(as.name("["),quote(sim$Data)), indTasForAux1))
                  callSim2 <- as.call(c(list(as.name("["),quote(sim$Data)), indTasForAux2))
                  indCalmWind <- which((eval(callSim1) == 0) & !is.na(sim$Data[indTasForAux]))
                  calmWind <- sim$Data[indTasForAux[indCalmWind,]]
                  sim$Data[indTasForAux] <- (eval(callSim2)*simTas$Data)/eval(callSim1)
                  if (length(calmWind)>0){
                        sim$Data[indTasForAux[indCalmWind,]] <- calmWind
                  }
            }
      }
      #       case {'prsn';'snowfall';'nieve'},
      #       if isempty(Pr),error('Precipitation is necessary for the correction of the snowfall');end
      #       prC=isimip(Pr.O,Pr.P,Pr.F,'variable','precipitation','datesobs',datesObs,'datesfor',datesFor,'threshold', threshold);
      #       indC=find(~isnan(Pr.F) & Pr.F>0);F(indC)=(F(indC).*prC(indC))./Pr.F(indC);
      if (!is.na(any(match(obs$Variable$varName,c("prsn","snowfall","nieve"))))) {
            indTas <- which(!is.na(match(obs$Variable$varName,c("pr","tp","precipitation","precip"))))
            if (length(indTas) == 0){
                  stop("Precipitation is needed to correct snow")
            }else{
                  if (!is.na(any(match(obs$Variable$varName,c("prsn","snowfall","nieve"))))) {
                        indTasMax <- which(!is.na(match(obs$Variable$varName,c("prsn","snowfall","nieve"))))
                  }
                  obsTas <- obs
                  indTasAux <- rep(list(bquote()), length(dimObs))
                  indTasAux[[obs.var.index]] <- indTas
                  callTas <- as.call(c(list(as.name("["),quote(obs$Data)), indTasAux))
                  obsTas$Data <- eval(callTas)
                  attr(obsTas$Data, "dimensions") <- attributes(obs$Data)$dimensions[setdiff(1:length(attributes(obs$Data)$dimensions),obs.var.index)]
                  obsTas$Variable$varName <- obs$Variable$varName[indTas]
                  obsTas$Dates$start <- obs$Dates[[indTas]]$start
                  obsTas$Dates$end <- obs$Dates[[indTas]]$end
                  prdTas <- pred
                  indTasAux <- rep(list(bquote()), length(dimPred))
                  indTasAux[[pred.var.index]] <- indTas
                  callTas <- as.call(c(list(as.name("["),quote(pred$Data)), indTasAux))
                  prdTas$Data <- eval(callTas)
                  attr(prdTas$Data, "dimensions") <- attributes(pred$Data)$dimensions[setdiff(1:length(attributes(pred$Data)$dimensions),pred.var.index)]
                  prdTas$Variable$varName <- pred$Variable$varName[indTas]
                  prdTas$Dates$start <- pred$Dates[[indTas]]$start
                  prdTas$Dates$end <- pred$Dates[[indTas]]$end
                  simTas <- sim
                  indTasAux <- rep(list(bquote()), length(dimFor))
                  indTasAux[[sim.var.index]] <- indTas
                  callTas <- as.call(c(list(as.name("["),quote(sim$Data)), indTasAux))
                  simTas$Data <- eval(callTas)
                  attr(simTas$Data, "dimensions") <- attributes(sim$Data)$dimensions[setdiff(1:length(attributes(sim$Data)$dimensions),sim.var.index)]
                  simTas$Variable$varName <- sim$Variable$varName[indTas]
                  simTas$Dates$start <- sim$Dates[[indTas]]$start
                  simTas$Dates$end <- sim$Dates[[indTas]]$end
                  simTas <- isimip(obsTas, prdTas, simTas, type = "multiplicative")# isimip
                  indTasFor <- rep(list(bquote()), length(dimFor))
                  for (d in 1:length(dimFor)){
                        indTasFor[[d]] <- 1:dimFor[d]
                  }
                  indTasForAux <- as.matrix(expand.grid(indTasFor))
                  indTasForAux1 <- rep(list(bquote()), length(dimFor))
                  indTasForAux2 <- rep(list(bquote()), length(dimFor))
                  indTasForAux1[[sim.var.index]] <- indTas
                  indTasForAux2[[sim.var.index]] <- indTasMax
                  callSim1 <- as.call(c(list(as.name("["),quote(sim$Data)), indTasForAux1))
                  callSim2 <- as.call(c(list(as.name("["),quote(sim$Data)), indTasForAux2))
                  indCalmWind <- which((eval(callSim1) == 0) & !is.na(sim$Data[indTasForAux]))
                  calmWind <- sim$Data[indTasForAux[indCalmWind,]]
                  sim$Data[indTasForAux] <- (eval(callSim2)*simTas$Data)/eval(callSim1)
                  if (length(calmWind)>0){
                        sim$Data[indTasForAux[indCalmWind,]] <- calmWind
                  }
            }
      }
      attr(sim$Data, "correction") <- "ISI-MIP"
      return(sim)
}
