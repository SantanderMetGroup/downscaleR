#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
loadGridDataset <- function(dataset, var, dictionary = NULL, lonLim = NULL, latLim = NULL, level = NULL, season = NULL, years = NULL) {
    if (is.null(dictionary) == FALSE) {
        shortName <- dictionaryLookup(dictionary, var)$ShortName
        standardName <- dictionaryLookup(dictionary, var)$StandardName
        isStandard <- TRUE
    } else {
        shortName <- var
        standardName <- var
        isStandard <- FALSE
    }
    # Open dataset
    gds <- J("ucar.nc2.dt.grid.GridDataset")$open(dataset)
    # Open grid
    grid <- gds$findGridByShortName(shortName)
    if (is.null(grid)) {
        stop("Variable requested not found. Check variable nomenclature")
    }
	# Get coordinate system
    gcs <- grid$getCoordinateSystem()
    # Retrieve geographical window / point location
    latLon <- getLatLonDomain(gcs, lonLim, latLim)
    # Retrieve time domain
    timeString <- unlist(strsplit(gsub("\\[|]|\\s", "", gcs$getTimeAxis()$getCalendarDates()$toString()), split=","))
    timeDates <- as.POSIXlt(timeString, format = "%Y-%m-%dT%H:%M:%SZ") ; rm(timeString)
    timePars <- getTimeDomain(timeDates, season, years)
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
    Data <- t(t(Data)[order(latLon$Grid[,2], latLon$Grid[,1]), ])
    gds$close()
    # Var transformation
    if (!is.null(dictionary)) {
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
#
#dataset = "/home/juaco/workspace/globalWildfires/fwi.ncml"
#var = "fwi"
#dictionary = NULL
#lonLim = c(-130,-80)
#latLim = c(20,50)
#level = NULL
#season = 6:9
#years = 1961:1965
#ex <- loadGridDataset(dataset, var, dictionary, lonLim, latLim, level, season, years)
#df <- cbind.data.frame(ex$LonLatCoords, colMeans(ex$Data))
# require(sp)
# require(maptools)
# data(wrld_simpl)
# l <- as(wrld_simpl, "SpatialLines")
# l1 <- list("sp.lines", l)
#coordinates(df) <- c(1,2)
#gridded(df) <- TRUE
#spplot(df, sp.layout = list(l1), col.regions = topo.colors(21), scales = list(draw = TRUE))

# Examples
# dataset = "/home/juaco/temp/meteoR/datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.ncml"
# var = "ta"
# dictionary = "/home/juaco/temp/meteoR/datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.dic"
# dataset <- "/home/juaco/temp/NCEP/NCEP_levels.ncml"
# di <- dataInventory(dataset)
# str(di)
# dictionary <- "/home/juaco/temp/NCEP/dictionary.csv"
# var = "uwnd"
# level = 850
# season = 5:7
# years = 2009
# lonLim = c(-10, 15)
# latLim = c(35,45)
# dictionary = NULL
####
#ex <- loadGridDataset(dataset = "/home/juaco/temp/NCEP/NCEP_levels.ncml", var = var, dictionary = dictionary, lonLim = c(-120, 1200), latLim = c(-10,75), level = 850, season = c(12,1,2), years = 2009:2010)
# ex$Dates$Start
# 
# 
# traceback()
# str(ex)
# require(sp)
# require(maptools)
# data(wrld_simpl)
# l <- as(wrld_simpl, "SpatialLines")
# l1 <- list("sp.lines", l)
# ex$LonLatCoords
# df <- cbind.data.frame(ex$LonLatCoords, colMeans(ex$Data))
# coordinates(df) <- c(1,2)
# gridded(df) <- TRUE
# spplot(df, sp.layout = list(l1), col.regions = topo.colors(21), scales = list(draw = TRUE))
