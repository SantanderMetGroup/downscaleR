# Examples
# NCEP
# dataset = "/home/juaco/temp/meteoR/datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.ncml"
# dictionary = "/home/juaco/temp/meteoR/datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.dic"
# predictors = c("ta", "zg", "hus")
# level = c(NULL, 850, 850)
# S4
dataset = "http://www.meteo.unican.es/tds5/dodsC/system4/System4_Seasonal_15Members.ncml"
dictionary = '/home/juaco/workspace/gitRepo/useRportal/dictionaries/System4_Seasonal_15Members.dic'
# dic <- read.csv(dictionary)
# dic
predictors = c("tas", "pr", "psl")
time.aggr = list("mean", "sum", 12)
level = NULL

di <- dataInventory(dataset)

str(di)



season = 5:7
years = 2009
lonLim = c(-10, 15)
latLim = c(35,45)




loadPredictors.GridDataset <- function(dataset, predictors, time.aggr, dictionary = NULL, lonLim = NULL, latLim = NULL, level = NULL, season = NULL, years = NULL) {
      dic <- read.csv(dictionary)
      
#       dic <- rbind.data.frame(dic, c("psl", "Mean_sea_level_pressure_surface", "24h", 0, 24, "mean", 0.00, 1, 0))

      if (is.null(level)) {
            level <- rep(NULL, length(predictors))
      }
      
      if (length(predictors) != length(time.aggr) | length(predictors) != length(time.aggr)) {
            stop("Differing lengths of 'predictors', 'levels' and/or 'time.aggr' arguments")
      }
      
      predictors.list <- list()
      
      for (p in 1:length(predictors)) {
     
            p=3
            
            var <- predictors[p]
            agg <- time.aggr[[p]]
      
                      
            if (is.character(agg)) {
                  dicRow <- which(dic$identifier == var & dic$aggr_fun == agg)
                  
            }
            if (is.numeric(agg)) {
                  dicRow <- which(dic$identifier == var & dic$aggr_fun == "none")
            }
            if (length(dicRow) == 0) {
                  stop("Requested time aggregation for predictor '", var, "' not available")
            }
            
            lev <- level[p]
            
            shortName <- as.character(dic$short_name[dicRow])
            standardName <- as.character(dic$identifier[dicRow])
            
            
            # Open dataset
            gds <- J("ucar.nc2.dt.grid.GridDataset")$open(dataset)
            # Open grid
            grid <- gds$findGridByShortName(shortName)
            if (is.null(grid)) {
                  stop("Variable requested not found. Check variable nomenclature")
            }
            # Get coordinate system
            gcs <- grid$getCoordinateSystem()
            
            a <- gcs$getTimeAxis()
            
            
            if (p == 1) {
                  # Retrieve geographical window / point location
                  latLon <- getLatLonDomain(gcs, lonLim, latLim)
                  # Retrieve time domain
                  timeString <- unlist(strsplit(gsub("\\[|]|\\s", "", gcs$getTimeAxis()$getCalendarDates()$toString()), split=","))
                  timeDates <- as.POSIXlt(timeString, format = "%Y-%m-%dT%H:%M:%SZ") ; rm(timeString)
                  timePars <- getTimeDomain(timeDates, season, years)
            }
            
            # Vertical levels
            vAxis <- tryCatch({gcs$getVerticalAxis()$getCoordValues()},
                        error = function(er) {
                              err <- NULL
                              return(err)
                        })
            if (is.null(vAxis) == FALSE) {
                  if (is.null(level)) {
                        stop("Variable with vertical levels: Argument 'level' required.\nSee data inventory for valid argument values")
                  }
                  levelInd <- match(level, vAxis) - 1
                  if (is.na(levelInd)) {
                        stop("Vertical level not found.\nCheck data inventory for valid vertical level values")
                  }
            }
      
            
            # Data retrieval
      message(paste("[",Sys.time(),"] - Retrieving data...", sep=""))
      dims <- unlist(strsplit(grid$getVariable()$getDimensionsString(), split = " "))
      Data <- matrix(ncol = nrow(latLon$Grid), nrow = length(unlist(timePars$timeIndList)))
      start <- rep(NA, length(dims))
      count <- start
      if (is.null(vAxis) == FALSE) {
            start[grep("^level", dims)] <- levelInd
            count[grep("^level", dims)] <- 1
      }
      start[grep("^lat", dims)] <- latLon$LatIndex[1]
      count[grep("^lat", dims)] <- length(latLon$LatIndex)
      indRowRange <- c(1)
      for (i in 1:length(timePars$timeIndList)) {
            DataRowRange <- seq.int(indRowRange, length(timePars$timeIndList[[i]]) + indRowRange - 1)
            indRowRange <- DataRowRange[length(DataRowRange)] + 1
            start[grep("^time", dims)] <- timePars$timeIndList[[i]][1]
            count[grep("^time", dims)] <- length(timePars$timeIndList[[i]])
            if (is.null(latLon$DatelineCrossInd)) {
                  start[grep("^lon", dims)] <- latLon$LonIndex[1]
                  count[grep("^lon", dims)] <- length(latLon$LonIndex)
                  aux <- grid$getVariable()$read(as.integer(start), as.integer(count))$getStorage()
                  aux.mat <- matrix(aux, ncol = nrow(latLon$Grid), byrow = TRUE) ; rm(aux)
                  Data[DataRowRange, ] <- aux.mat ; rm(aux.mat)
            } else {
                  lonStartList <- list(1 : (latLon$DatelineCrossInd), (latLon$DatelineCrossInd + 1) : length(latLon$LonIndex))
                  DataColRangeList <- list(1 : (latLon$DatelineCrossInd * length(latLon$LatIndex)), (latLon$DatelineCrossInd * length(latLon$LatIndex) + 1) : nrow(latLon$Grid))
                  for (j in 1:length(lonStartList)) {
                        start[grep("^lon", dims)] <- latLon$LonIndex[lonStartList[[j]][1]]
                        count[grep("^lon", dims)] <- length(lonStartList[[j]])
                        aux <- grid$getVariable()$read(as.integer(start), as.integer(count))$getStorage()
                        aux.mat <- matrix(aux, ncol = length(DataColRangeList[[j]]), byrow = TRUE) ; rm(aux)
                        Data[DataRowRange, DataColRangeList[[j]]] <- aux.mat ; rm(aux.mat)
                  }
                  latLon$Grid <- rbind(latLon$Grid[which(latLon$Grid[ ,1] < 0), ], latLon$Grid[which(latLon$Grid[ ,1] >= 0), ])
            }
      }
      gds$close()
      # Var transformation
      if (is.null(dictionary) == FALSE) {
            l <- dictionaryTransform(Data, shortName, dictionary, timePars)
            for (i in 1:length(l)) {
                  assign(names(l)[i], l[[i]])
            } ; rm(l)
      } else {
            varTimeStep <- difftime(timePars$dateSlice[2], timePars$dateSlice[1])
            dateSliceStart <- timePars$dateSlice
            dateSliceEnd <- timePars$dateSlice + varTimeStep
      }
      dateSliceStart <- as.POSIXlt(dateSliceStart)
      dateSliceEnd <- as.POSIXlt(dateSliceEnd)
      message(paste("[",Sys.time(),"] - Done.", sep=""))
      return(list("VarName" = standardName, "isStandard" = isStandard, "Level" = level, "Dates" = list("Start" = dateSliceStart,"End" = dateSliceEnd), "LonLatCoords" = latLon$Grid, "Data" = Data))
}
# End
