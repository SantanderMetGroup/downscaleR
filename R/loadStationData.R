#' @title Load station data
#' @description Load observations data from station datasets in standard ASCII format.

#' @param source.dir A valid path to the directory containing the station files
#' @param file.format Wether the stations data are stored in a netCDF or ASCII (default) file. See details for standard format definition.
#' @param var Character string indicating the variable to be loaded. Note that the notation depends on wether the dictionary is being used or not.
#' @param dictionary Optional. A full path to the file containing the dictionary. See details
#' @param stationID Optional. A character vector indicating the code names of the stations to be loaded.
#' @param lonLim Vector of length = 2, with minimum and maximum longitude coordinates of the bounding box
#'  selected. For single-point queries, a numeric value with the longitude coordinate. See details.
#' @param latLim Same as \code{lonLim}, but for the selection of the latitudinal range. See details.
#' @param season An integer vector specifying the desired season (in months, January = 1 ..., December = 12).
#'  Options include one to several (contiguous) months. If NULL (default), full year selections is performed (same as \code{season = 1:12})
#' @param years Optional vector of years to select. Default (NULL) to all available years
#' @param tz A time zone specification to be used for the conversion of dates, if one is required
#' (i.e., if the time zone of the dataset does not correspond to the system-specific one; see
#' \code{\link[base]{timezones}} for details). The default assumes "GMT" (UTC, Universal Time, Coordinated),
#' but be very careful as station datasets are quite heterogeneous.
#' 
#' @return a list with the following elements:
#' \itemize{
#' \item \code{Variable}. Name of the variable
#' \item \code{Data}. Dates are ordered by rows and Stations in columns. Names are station codes
#' \item \code{xyCoords}. A 2-D matrix with longitude and latitudes of the stations
#' \item \code{Dates}. A list with the verification time interval of each record in the time series.
#'  This is represented by a list with two elements: \code{start} and \code{end}, representing the
#'  lower and upper bounds of the verification period
#' \item \code{Metadata}. A list of variable length depending on the available metadata associated
#' to each observation. If no metadata are provided, at least the station codes (compulsory) are displayed.
#' }
#' 
#' @references \url{https://github.com/SantanderMetGroup/downscaleR/wiki/Observation-Data-format} 
#' 
#' @export
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' 
#' @seealso \code{\link{loadGridData}}, \code{\link{dataInventory}}
#' @aliases loading

loadStationData <- function(source.dir, file.format = c("ascii", "netcdf"), var, 
            stationID = NULL, lonLim = NULL, latLim = NULL, season = NULL,
            years = NULL, tz = "GMT") {
            file.format <- match.arg(file.format, choices = c("ascii", "netcdf"))
      if ((!is.null(lonLim) | !is.null(latLim)) & !is.null(stationID)) { 
            lonLim <- NULL 
            latLim <- NULL
            warning("lonLim/latLim arguments ignored as Station Codes have been specified.")
      }
      if (file.format == "ascii") {
            out <- loadStationData.ASCII(source.dir, var, stationID, lonLim, latLim, season, years)
      } else if (file.format == "netcdf") {
            out <- loadStationData.NetCDF(source.dir, var, stationID, lonLim, latLim, season, years)
      }
      varTimeStep <- difftime(out$Dates[2], out$Dates[1])
      dateSliceStart <- as.POSIXct(out$Dates)
      dateSliceEnd <- as.POSIXct(as.POSIXlt(out$Dates + varTimeStep))
      out$Dates <- list("start" = dateSliceStart, "end" = dateSliceEnd)
      message(paste("[", Sys.time(), "] Done.", sep = ""))
      return(out)
}
# End      





