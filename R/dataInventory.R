#' @title Dataset inventory
#' @description Function to provide a quick overview of a climate dataset
#'  (either stations or gridded data)
#' @param dataset A character string poiting to the target. Either a directory containing the data
#'  in the case of station data in standard ASCII format (see \code{\link{loadStationData}}) for details),
#'  or a target file (a NcML) in the case of other types of gridded data (reanalysis, gridded observations ...,
#'  see \code{\link{loadGridData}} for details).
#' @param return.stats Optional logical flag indicating if summary statistics of the dataset
#'  should be returned with the inventory. Only used for station data.
#'  
#' @return A list of components describing the variables and other characteristics of the target dataset.
#' 
#' @note The variable names returned correspond to the original names of the variables as stored in the dataset,
#' and not to the standard naming convention defined in the vocabulary.
#' 
#' @examples \donttest{
#' gsn <- file.path(find.package("downscaleR"), "datasets/observations/GSN_Iberia")
#' di <- dataInventory(gsn)
#' str(di)
#' # To obtain summary statistics of the variables stored:
#' di.stats <- dataInventory(gsn, return.stats = TRUE)
#' print(di.stats$Summary.stats)
#' } 
#' 
#' @seealso \code{\link{stationInfo}} for a quick overview of available stations in station datasets.
#' @export
#' @author J Bedia \email{joaquin.bedia@@gmail.com}

dataInventory <- function(dataset, return.stats = FALSE) {
      rs <- return.stats
      message(paste("[", Sys.time(), "] Doing inventory ...", sep = ""))
      if (isTRUE(file.info(dataset)$isdir)) {
            out <- dataInventory.ASCII(dataset, rs)
      } else {
            out <- dataInventory.NetCDF(dataset)
      }
      message(paste("[", Sys.time(), "] Done.", sep = ""))
      return(out)
}
# End

#' @title Retrieve station info
#' @description Get a quick overview of the stations contained in a stations dataset
#' @param dataset A character string indicating the database to be accessed. For station data in standard ASCII format,
#' this is the path to the directory the dataset lives in.
#' @param plot Logical indicating if a map displaying the station locations should be displayed. Default to TRUE.
#' @return A data.frame with the station codes in its first column, followed by longitude and latitude,
#'  and the rest of metadata (if any) in the following columns. If \code{plot = TRUE}, also a map of the station locations.
#' @note For an adequate map display the station coordinates must be in decimal degrees.
#' @seealso \code{\link{dataInventory}} to obtain a more exhaustive report the dataset.
#' @export
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @examples \donttest{
#' gsn <- file.path(find.package("downscaleR"), "datasets/observations/GSN_Iberia")
#' print(stationInfo(gsn))
#' } 
#' 

stationInfo <- function(dataset, plot = TRUE) {
      di <- dataInventory(dataset, return.stats = FALSE)
      ll <- di$Stations$LonLatCoords
      l <- lapply(1:length(di$Stations$other.metadata), function(x) {di$Stations$other.metadata[[x]]})
      df <- do.call("cbind.data.frame", l)
      l <- NULL
      names(df) <- names(di$Stations$other.metadata)
      stids <- di$Stations$station_id
      df <- cbind.data.frame("stationID" = stids,"longitude" = ll[ ,1], "latitude" = ll[ ,2], df)
      rownames(df) <- NULL
      if (isTRUE(plot)) {
            x.off <- diff(range(ll[ ,1]))*.3
            y.off <- diff(range(ll[ ,2]))*.3
            x.ran <- c(min(ll[ ,1]) - x.off, max(ll[ ,1]) + x.off)
            y.ran <- c(min(ll[ ,2]) - y.off, max(ll[ ,2]) + y.off)
            plot(ll, asp = 1, xlim = x.ran, ylim = y.ran, col = "blue", pch = 10)
            world(add = TRUE)
            text(x = ll[,1], y = ll[,2], labels = stids, pos = 3, cex = .7, col = "red")
      }
      return(df)
}
# End



