#' Load observations data

#' Loads observations data from station and gridded datasets, either in standard ASCII format or netCDF

#' @param source.dir A valid path to the directory containing the station files
#' @param file.format Wether the stations data are stored in a netCDF or ASCII (default) file. See details for standard format definition.
#' @param var Character string indicating the variable to be loaded. Note that the notation depends on wether the dictionary is being used or not.
#' @param dictionary Optional. A full path to the file containing the dictionary. See details
#' @param stationID Optional. A character vector indicating the code names of the stations to be loaded.
#' @param lonLim
#' @param latLim
#' @param season
#' @param years
#' @return a list with the following elements:
#' \begin{itemize}
#' \item{VarName}{name of the variable requested (Note that the name can be either standard or not depending on wether a dictionary is being used)}
#' Stations A data frame containing the station codes and their names (locations).
#' LonLatCoords: A 2-D matrix with longitude and latitudes of the stations
#' Dates: A list with the verification times of each record in the time series. This is represented by a list with two elements:
	# 'Start': The lower time bound of the verification period
	# 'End': The upper time bound of the verification period
#' Data A matrix with the values of the variable. Dates are ordered by rows and Stations in columns.
#' \end{itemize}
#' @export
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @seealso \code{\link{loadGridDataset}}, \code{\link{dataInventory}}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# source.dir = "/home/DATA/CLIMATE/GSN_World/"
# file.format = "ascii"
# var = "tasmax"
# dictionary = "/home/DATA/CLIMATE/GSN_World/GSN_World.csv"
# stationID = NULL
# lonLim = c(-5,10)
# latLim = c(37,43)
# season = 6:8
# years = 1990:1995

loadObservations <- function(source.dir, file.format = c("ascii", "netcdf"), var, stationID = NULL, lonLim = NULL, latLim = NULL, season = NULL, years = NULL) {
    file.format <- match.arg(file.format, choices = c("ascii", "netcdf"))
    if ((!is.null(lonLim) | !is.null(latLim)) & !is.null(stationID)) { 
        lonLim <- NULL 
        latLim <- NULL
        warning("lonLim/latLim arguments ignored as Station Codes have been specified.")
    }
    if (file.format == "ascii") {
        out <- loadObservations.ASCII(source.dir, var, stationID, lonLim, latLim, season, years)
    } else if (file.format == "netcdf") {
        out <- loadObservations.NetCDF(source.dir, var, stationID, lonLim, latLim, season, years)
    }
    varTimeStep <- difftime(out$time[2], out$time[1])
    dateSliceStart <- out$time
    dateSliceEnd <- as.POSIXlt(out$time + varTimeStep)
    out$time <- list("Start" = dateSliceStart, "End" = dateSliceEnd)
    message(paste("[", Sys.time(), "] Done.", sep = ""))
    return(out)
}
# End      





