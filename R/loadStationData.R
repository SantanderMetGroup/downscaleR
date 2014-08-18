#' @title Load station data
#' @description Load observations data from station datasets in standard ASCII format.
#'
#' @template templateParams 
#' @param file.format Wether the stations data are stored in a netCDF or ASCII (default) file. 
#' Currently only the standard \dQuote{ASCII} format supported. See references for details.
#' @param stationID Optional. A character vector indicating the code names of the stations to be loaded.
#' @param tz A time zone specification to be used for the conversion of dates, if one is required
#' (i.e., if the time zone of the dataset does not correspond to the system-specific one; see
#' \code{\link[base]{timezones}} for details). The default assumes "GMT" (UTC, Universal Time, Coordinated),
#' but be very careful as station datasets are quite heterogeneous.
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
            years = NULL, tz = "GMT") {
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
      out$Dates <- list("start" = dateSliceStart, "end" = dateSliceEnd)
      attr(out$Data, "dimensions") <- c("time", "station")
      message(paste("[", Sys.time(), "] Done.", sep = ""))
      return(out)
}
# End      





