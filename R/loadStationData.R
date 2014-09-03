#' @title Load station data
#' @description Load observations data from station datasets in standard ASCII format.
#'
#' @template templateParams 
#' @param file.format Wether the stations data are stored in a netCDF or ASCII (default) file. 
#' Currently only the standard \dQuote{ASCII} format supported. See references for details.
#' @param stationID Optional. A character vector indicating the code names of the stations to be loaded.
#' @param tz A time zone specification to be used for the conversion of dates, if one is required
#' (i.e., if the time zone of the dataset does not correspond to the system-specific one; see
#' \code{\link[base]{timezones}} for details). Default to unspecified (i.e. \code{tz = ""}).
#' 
#' @return a list with the following elements:
#' \itemize{
#' \item \code{Variable}. Name of the variable
#' \item \code{Data}. A 2-D matrix containing the data. Dates are ordered by rows and Stations by columns, 
#' following the order indicated in the \code{Metadata}.
#' \item \code{xyCoords}. A 2-D matrix with longitude and latitudes of the stations
#' \item \code{Dates}. A list with the verification time interval of each record in the time series.
#'  This is represented by a list with two elements: \code{start} and \code{end}, representing the
#'  lower and upper bounds of the verification period
#' \item \code{Metadata}. A list of variable length depending on the available metadata associated
#' to each observation. If no metadata are provided, at least the station codes (compulsory) are displayed.
#' }
#' 
#' @template templateGeolocation
#' 
#' @note Unlike gridded datasets, station data do not use a dictionary for variable homogenization. Thus, users
#' must take care of variable units and eventual conversions when necessary.
#' 
#' @references \url{https://github.com/SantanderMetGroup/downscaleR/wiki/Observation-Data-format} 
#' 
#' @export
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' 
#' @family loading

loadStationData <- function(dataset, file.format = "ascii", var, 
            stationID = NULL, lonLim = NULL, latLim = NULL, season = NULL,
            years = NULL, tz = "") {
            file.format <- match.arg(file.format, choices = c("ascii", "netcdf"))
      if ((!is.null(lonLim) | !is.null(latLim)) & !is.null(stationID)) { 
            lonLim <- NULL 
            latLim <- NULL
            warning("lonLim/latLim arguments ignored as Station Codes have been specified.")
      }
      if (file.format == "ascii") {
            out <- loadStationData.ASCII(dataset, var, stationID, lonLim, latLim, season, years, tz)
      } else {
            stop("Unrecognized/not supported station data format")
      }
#       } else if (file.format == "netcdf") {
#             out <- loadStationData.NetCDF(dataset, var, stationID, lonLim, latLim, season, years, tz)
#       }
      varTimeStep <- difftime(out$Dates[2], out$Dates[1])
      dateSliceStart <- as.POSIXct(out$Dates)
      dateSliceEnd <- as.POSIXct(as.POSIXlt(out$Dates + varTimeStep))
      usetz <- ifelse(identical(tz, ""), FALSE, TRUE)
      dateSliceStart <- format.POSIXct(dateSliceStart, "%Y-%m-%d %H:%M:%S", usetz = usetz)
      dateSliceEnd <- format.POSIXct(dateSliceEnd, "%Y-%m-%d %H:%M:%S", usetz = usetz)
      out$Dates <- list("start" = dateSliceStart, "end" = dateSliceEnd)
      attr(out$Data, "dimensions") <- c("time", "station")
      message(paste("[", Sys.time(), "] Done.", sep = ""))
      return(out)
}
# End      
################################################################################


#' Load station data in standard ASCII format
#' 
#' Load station data in standard ASCII format, being the standard defined in the framework
#' of the action COST VALUE
#' 
#' @param dataset Directory containing the stations dataset
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

loadStationData.ASCII <- function(dataset, var, stationID, lonLim, latLim, season, years, tz) {
      aux <- read.csv(file.path(dataset, "stations.txt"), stringsAsFactors = FALSE, strip.white = TRUE)
      # Station codes
      stids <- read.csv(file.path(dataset, "stations.txt"), colClasses = "character")[ ,grep("station_id", names(aux), ignore.case = TRUE)]
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
      fileInd <- grep(paste("^", var, "\\.txt", sep = ""), list.files(dataset))
      if(length(fileInd) == 0) {
            stop("[", Sys.time(),"] Variable requested not found")
      }
      timeString <- read.csv(list.files(dataset, full.names = TRUE)[fileInd], colClasses = "character")[ ,1]
      if (nchar(timeString[1]) == 8) {
            timeDates <- strptime(timeString, "%Y%m%d", tz = tz)  
      }
      if (nchar(timeString[1]) == 10) {
            timeDates <- strptime(timeString, "%Y%m%d%H", tz = tz)
      }
      timeString <- NULL
      timePars <- getTimeDomainStations(timeDates, season, years)
      ## missing data code
      vars <- read.csv(list.files(dataset, full.names=TRUE)[grep("variables", list.files(dataset), ignore.case = TRUE)])
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
      Data <- unname(as.matrix(read.csv(list.files(dataset, full.names = TRUE)[fileInd], na.strings = na.string)[timePars$timeInd, stInd + 1]))
      #       names(Data) <- stids
      ## Metadata
      message("[", Sys.time(), "] Retrieving metadata ...", sep = "")
      # Assumes that at least station ids must exist, and therefore meta.list is never empty
      ind.meta <- c(1:length(names(aux)))[-pmatch(c("longitude", "latitude"), names(aux))]
      meta.list <- as.list(aux[stInd,ind.meta])
      aux <- NULL  
      return(list("Variable" = list("varName" = var), "Data" = Data, "xyCoords" = coords, "Dates" = timePars$timeDates, "Metadata" = meta.list))
}
# End




