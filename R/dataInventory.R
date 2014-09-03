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
################################################################################

#' @title Data inventory of standard ASCII station datasets
#' @description Make data inventory from a station dataset in standard ASCII format
#' @param dataset path to the directory containng the dataset
#' @param rs Logical. return stats?
#' @return A data inventory 
#' @author J Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal

dataInventory.ASCII <- function(dataset, rs) {
      lf <- list.files(dataset, full.names = TRUE)
      stations <- read.csv(lf[grep("stations", lf, ignore.case = TRUE)], strip.white = TRUE, stringsAsFactors = FALSE)
      vars <- read.csv(lf[grep("variables", lf, ignore.case = TRUE)], strip.white = TRUE, colClasses = "character")
      # Var info 
      var.info <- list("variable" = vars[ ,grep("variable", names(vars), ignore.case = TRUE)], "longname" = vars[ ,grep("longname", names(vars), ignore.case = TRUE)], "unit" = vars[ ,grep("unit", names(vars), ignore.case = TRUE)], "missing.code" = vars[ ,grep("missing_code", names(vars), ignore.case = TRUE)])
      var.info <- do.call("cbind.data.frame", var.info)
      # Station info
      timeString <- read.csv(lf[grep(paste("^", vars[ ,1][1], "\\.*", sep = ""), list.files(dataset))], colClasses = "character")[ ,1]
      if (nchar(timeString[1]) == 8) {
            timeDates <- strptime(timeString, "%Y%m%d")  
      } else {
            timeDates <- strptime(timeString, "%Y%m%d%H")
      }
      rm(timeString)
      timeAxis <- list("startDate" = min(timeDates), "endDate" = max(timeDates), "timeStep" = difftime(timeDates[2], timeDates[1], units = "h"))    
      rm(timeDates)
      station_id <- as.character(stations[ ,grep("station_id", names(stations), ignore.case = TRUE)])
      lon <- stations[ ,grep("^longitude$", names(stations), ignore.case = TRUE)]
      lat <- stations[ ,grep("^latitude$", names(stations), ignore.case = TRUE)]
      LonLatCoords <- cbind(lon, lat)
      rownames(LonLatCoords) <- station_id
      rm(lon, lat)
      other.metadata <- as.list(stations[ ,-pmatch(c("station_id", "longitude", "latitude"), names(stations))])
      station.info <- list("station_id" = station_id, "LonLatCoords" = LonLatCoords, "times" = timeAxis, "other.metadata" = other.metadata)
      rm(station_id, timeAxis, LonLatCoords, other.metadata)
      info <- list("Stations" = station.info, "Variables" = var.info, "Summary.stats" = NULL)
      if (isTRUE(rs)) {
            na.perc <- function(x) {round(length(which(is.na(x))) * 100 / length(x), 1)}
            aux.mat <- matrix(ncol = nrow(vars), nrow = length(station.info$station_id), dimnames = list(station.info$station_id, vars[ ,grep("variable", names(vars), ignore.case = TRUE)]))
            stats.list <- list("missing.percent" = aux.mat, "min" = aux.mat, "max" = aux.mat, "mean" = aux.mat)
            for (i in 1:nrow(vars)) {
                  var <- read.csv(lf[grep(paste("^", var.info[i,1], "\\.*", sep = ""), list.files(dataset))], na.strings = var.info$missing.code[i])[ ,-1]
                  stats.list[[1]][ ,i] <- apply(var, 2, na.perc)
                  stats.list[[2]][ ,i] <- apply(var, 2, min, na.rm = TRUE)
                  stats.list[[3]][ ,i] <- apply(var, 2, max, na.rm = TRUE)
                  stats.list[[4]][ ,i] <- apply(var, 2, mean, na.rm = TRUE)
                  rm(var)
            }
            info <- list("Stations" = station.info, "Variables" = var.info, "Summary.stats" = stats.list)
      }
      return(info)
}
# End
################################################################################

#' @title Inventory of a gridded dataset
#' 
#' @description Returns a list with summary information about the variables stored in a gridded dataset.
#' Sub-routine of \code{dataInventory}
#' 
#' @param dataset A full path to the file describing the dataset (NcML)
#' 
#' @return A (named) list whose length is determined by the number of variables stored in the dataset,
#' its names corresponding to the short names of the variables.
#' For each variable, information on the variable long name, data type, units and
#' characteristics of its dimensions is provided.
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' 
#' @keywords internal
#' 
#' @import rJava



dataInventory.NetCDF <- function(dataset) {
      gds <- J("ucar.nc2.dt.grid.GridDataset")$open(dataset)
      varNames <- unlist(strsplit(gsub("\\[|]|\\s", "", gds$getGrids()$toString()), ","))
      rm.ind <- grep("^lon|^lat", varNames)  
      if (length(rm.ind) > 0) {
            varNames <- varNames[-rm.ind]
      }
      if (length(varNames) == 0) {
            stop("No variables found")
      } else {
            var.list <- list()
            for (i in 1:length(varNames)) {
                  message("[", Sys.time(), "] Retrieving info for \'", varNames[i], "\' (", length(varNames) - i, " vars remaining)")
                  description <- gds$getDataVariable(varNames[i])$getDescription()
                  varName <- gds$getDataVariable(varNames[i])$getShortName()
                  dataType <- gds$getDataVariable(varNames[i])$getDataType()$toString()
                  units <- gds$getDataVariable(varNames[i])$getUnitsString()
                  grid <- gds$findGridByShortName(varName)
                  dim.list <- scanVarDimensions(grid)
                  var.list[[i]] <- list("Description" = description, "DataType" = dataType, "Units" = units, "Dimensions" = dim.list)
            }
            names(var.list) <- varNames
      }
      gds$close()
      return(var.list)
}
# End
################################################################################


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



