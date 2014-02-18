loadObservations.ASCII <- function(source.dir, var, stationID = NULL, lonLim = NULL, latLim = NULL, season = NULL, years = NULL) {
      aux <- read.csv(paste(source.dir, "/stations.txt", sep = ""), stringsAsFactors = FALSE, strip.white = TRUE)
      # Station codes
      stids <- read.csv(paste(source.dir, "/stations.txt", sep = ""), colClasses = "character")[ ,grep("station_id", names(aux), ignore.case = TRUE)]
      if (!is.null(stationID)) {
		stInd <- match(stationID, stids)
		if (any(is.na(stInd))) {
			stop("'stationID' values not found.\nCheck data inventory")
		}
	} else {
		stInd <- 1:length(stids)
	}
      ## Longitude and latitude
	lons <- aux[ ,grep("^longitude$", names(aux), ignore.case = TRUE)]
	lats <- aux[ ,grep("^latitude$", names(aux), ignore.case = TRUE)]
	if (!is.null(lonLim)) {
		latLon <- getLatLonDomainStations(lonLim, latLim, lons, lats)
		if (length(latLon$stInd) == 0) {
			stop("No stations were found in the selected spatial domain")
		}
		stInd <- latLon$stInd
		coords <- latLon$stCoords
		rm(latLon)
	} else {
		coords <- matrix(cbind(lons, lats)[stInd, ], ncol = 2)
	}
      stids <- stids[stInd]
      dimnames(coords) <- list(stids, c("longitude", "latitude"))
      ## Time dimension
	fileInd <- grep(paste("^", var, "\\.txt", sep = ""), list.files(source.dir))
	if(length(fileInd) == 0) {
            stop("[", Sys.time(),"] Variable requested not found")
	}
      timeString <- read.csv(list.files(source.dir, full.names = TRUE)[fileInd], colClasses = "character")[ ,1]
      if (nchar(timeString[1]) == 8) {
		timeDates <- strptime(timeString, "%Y%m%d")  
	} else {
		timeDates <- strptime(timeString, "%Y%m%d%H")
	}
      rm(timeString)
      timePars <- getTimeDomain(timeDates, season, years)
      rm(timeDates)
      ## missing data code
      vars <- read.csv(list.files(source.dir, full.names=TRUE)[grep("variables", list.files(source.dir), ignore.case = TRUE)])
      miss.col <- grep("missing_code", names(vars), ignore.case = TRUE)
      if(length(miss.col) > 0) {
            na.string <- vars[grep(var, vars[ ,grep("variable", names(vars), ignore.case = TRUE)]), miss.col]
            rm(vars, miss.col)
      } else {
            na.string <- NA
      }
      # Data retrieval
      message(paste("[", Sys.time(), "] Loading data ...", sep = ""))
      Data <- as.data.frame(read.csv(list.files(source.dir, full.names = TRUE)[fileInd], na.strings = na.string)[unlist(timePars$timeIndList), stInd + 1])
      names(Data) <- stids
	## Metadata
      message(paste("[", Sys.time(), "] Retrieving metadata ...", sep = ""))
      ind.meta <- c(1:length(names(aux)))[-pmatch(c("station_id", "longitude", "latitude"), names(aux))]
      if (length(ind.meta) == 0) {
            meta.list <- NULL
      } else {
            meta.list <- as.list(aux[stInd,ind.meta])
      }
      rm(aux)  
      return(list("variable" = var, "station_id" = stids, "LonLatCoords" = SpatialPoints(coords), "time" = timePars, "metadata" = meta.list, "Data" = Data))
}
# End