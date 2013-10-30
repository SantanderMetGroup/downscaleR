# TODO: Add comment
# 
# Author: juaco
###############################################################################
loadObservations.ASCII <- function(source.dir, shortName, stationID = NULL, lonLim = NULL, latLim = NULL, season = NULL, years = NULL) {
	source.dir <- source.dir
	shortName <- shortName
	stationID <- stationID 
	lonLim <- lonLim
	latLim <- latLim
	season <- season
	years <- years
	## Station Ids
	stids <- read.csv(paste(source.dir, "/stations.txt", sep = ""), colClasses = "character")[ ,1]
	if (is.null(stationID) == FALSE) {
		stInd <- match(stationID, stids)
		if (any(is.na(stInd))) {
			stop("'stationID' values not found.\nCheck data inventory")
		}
	} else {
		stInd <- 1:length(stids)
	}
	## Longitude and latitude
	lons <- read.csv(paste(source.dir, "/stations.txt", sep = ""))[ ,3]
	lats <- read.csv(paste(source.dir, "/stations.txt", sep = ""))[ ,4]
	if (is.null(lonLim) == FALSE) {
		latLon <- getLatLonDomainStations(lonLim, latLim, lons, lats)
		if (length(latLon$stInd) == 0) {
			stop("No stations were found in the selected spatial domain")
		}
		stInd <- latLon$stInd
		coords <- latLon$stCoords
		rm(latLon)
	} else {
		coords <- cbind(lons, lats)[stInd, ]
	}
	if (is.null(nrow(coords)) == FALSE) {
		row.names(coords) <- stids[stInd]
	}
	## Time dimension
	fileInd <- grep(paste("^", shortName, "\\.txt", sep = ""), list.files(source.dir))
	timeString <- read.csv(list.files(source.dir, full.names=TRUE)[fileInd], colClasses = "character")[ ,1]
	if (nchar(timeString[1]) == 8) {
		timeDates <- strptime(timeString, "%Y%m%d")  
	} else {
		timeDates <- strptime(timeString, "%Y%m%d%H")
	}
	timePars <- getTimeDomain(timeDates, season, years)
	# Data retrieval
	Data <- as.data.frame(read.csv(list.files(source.dir, full.names = TRUE)[fileInd])[unlist(timePars$timeIndList), stInd + 1])
	names(Data) <- stids[stInd]
	stNames <- gsub("^\\s|\\s*$", "", read.csv(paste(source.dir, "/stations.txt", sep = ""))[stInd,2])
	stations <- cbind.data.frame("stationIDs" = stids[stInd], "stationNames" = stNames, stringsAsFactors = FALSE)
	return(list("VarName" = shortName, "Stations" = stations, "LonLatCoords" = coords, "TimePars" = timePars, "Data" = Data))
}
# End

# Example
#list.files('/home/juaco/temp/VALUE_validationTest_v2')
#source.dir = '/home/juaco/temp/VALUE_validationTest_v2'
#dictionary = "/home/juaco/temp/VALUE_validationTest_v2/dictionary.dic"
#var = "ta"
#lonLim = c(-10,10)
#latLim = c(35,44)
#ex <- loadObservations.ASCII(source.dir, shortName = "temp", stationID = NULL, , season = c(12,1,2), years = 1987:2000)
