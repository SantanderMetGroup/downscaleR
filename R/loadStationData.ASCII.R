#' Load station data in standard ASCII format
#' 
#' Load station data in standard ASCII format, being the standard defined in the framework
#' of the action COST VALUE
#' 
#' @param source.dir Directory containing the stations dataset
#' @param var Character string. Name of the variable, as defined in the dataset
#' @param stationID Character string. Optional, Id code(s) of the station(s) selected
#' @param lonLim numeric vector of length 1 or two defining X coordinate(s).
#'  See \code{\link{loadStationData}} for details.
#' @param latLim numeric vector of length 1 or two defining Y coordinate(s).
#'  See \code{\link{loadStationData}} for details.
#' @param season Numeric vector of months defining the season
#' @param years Numeric vector of years
#' @return A list with several components. See \code{\link{loadStationData}} for details.
#' @references \url{https://github.com/SantanderMetGroup/downscaleR/wiki/Observation-Data-format} 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal

loadStationData.ASCII <- function(source.dir, var, stationID = NULL, lonLim = NULL, latLim = NULL, season = NULL, years = NULL) {
      aux <- read.csv(file.path(source.dir, "stations.txt"), stringsAsFactors = FALSE, strip.white = TRUE)
      # Station codes
      stids <- read.csv(file.path(source.dir, "stations.txt"), colClasses = "character")[ ,grep("station_id", names(aux), ignore.case = TRUE)]
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
		latLon <- NULL
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
	      timeDates <- strptime(timeString, "%Y%m%d", tz = tz)  
	}
      if (nchar(timeString[1]) == 10) {
		timeDates <- strptime(timeString, "%Y%m%d%H", tz = tz)
	}
      timeString <- NULL
      timePars <- getTimeDomainStations(timeDates, season, years)
      ## missing data code
      vars <- read.csv(list.files(source.dir, full.names=TRUE)[grep("variables", list.files(source.dir), ignore.case = TRUE)])
      miss.col <- grep("missing_code", names(vars), ignore.case = TRUE)
      if(length(miss.col) > 0) {
            na.string <- vars[grep(var, vars[ ,grep("variable", names(vars), ignore.case = TRUE)]), miss.col]
            vars <- NULL
            miss.col <- NULL
      } else {
            na.string <- NA
      }
      # Data retrieval
      message("[", Sys.time(), "] Loading data ...", sep = "")
      Data <- as.data.frame(read.csv(list.files(source.dir, full.names = TRUE)[fileInd], na.strings = na.string)[timePars$timeInd, stInd + 1])
      names(Data) <- stids
	## Metadata
      message("[", Sys.time(), "] Retrieving metadata ...", sep = "")
      # Assumes that at least station ids must exist, and therefore meta.list is never empty
      ind.meta <- c(1:length(names(aux)))[-pmatch(c("longitude", "latitude"), names(aux))]
      meta.list <- as.list(aux[stInd,ind.meta])
      aux <- NULL  
      return(list("Variable" = var, "Data" = Data, "xyCoords" = coords, "Dates" = timePars$timeDates, "Metadata" = meta.list))
}
# End

