#' @title Retrieve dimension information for a gridded variable

#' @description Goes through variable dimensions and retrieves information for the inventory.
#' This is a sub-routine of \code{dataInventory.NetCDF}, \code{makeSubset}.
#' 
#' @import rJava
#' @importFrom downscaleR.java javaCalendarDate2rPOSIXlt
#' @importFrom downscaleR.java javaString2rChar
#' 
#' @param grid A java object of the class \sQuote{ucar.nc2.dt.grid.GeoGrid}
#' @return a list of length \emph{N}, being N the number of dimensions defining the grid shape.
#' @details All possible dimensions of a gridded datasets are searched and included if existing. This is done
#'  following the canonical dimension ordering [runtime, member, time, level, lat, lon], although function
#'  does not assume this ordering and therefore non-standard orderings should be adequately treated
#'  (this might be the case for instance when creating a new aggregation dimension via NcML).
#'  In case of existing vertical levels for a given dataset, all possible dataset level values are scanned
#'  one by one and those that do not exist for that particular variable are removed. This way, the inventory returns
#'  only the non-empty level values for each variable. In case of horizontal 2D axes,
#'  coordinate values are truncated to 1D for conciseness in the inventory.
#'  
#' @references \url{http://www.unidata.ucar.edu/software/thredds/v4.3/netcdf-java/v4.3/javadocAll/ucar/nc2/dt/grid/GeoGrid.html}
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @export
#' @keywords internal

scanVarDimensions <- function(grid) {
      gridShape <- grid$getShape()
      gcs <- grid$getCoordinateSystem()
      dim.list <- list()
      length(dim.list) <- length(gridShape)
      if (grid$getRunTimeDimensionIndex() >= 0) {
            rtDimIndex <- grid$getRunTimeDimensionIndex() + 1
            axis <- gcs$getRunTimeAxis()
            dim.list[[rtDimIndex]] <- list("Type" = axis$getAxisType()$toString(), "Units" = axis$getUnitsString(), "Values" = axis$getCoordValues())
            names(dim.list)[rtDimIndex] <- "runtime" # axis$getDimensionsString()
      }
      if (grid$getEnsembleDimensionIndex() >= 0) {
            ensDimIndex <- grid$getEnsembleDimensionIndex() + 1
            axis <- gcs$getEnsembleAxis()
            ens.values <- tryCatch(axis$getCoordValues(), error = function(er) {er <- javaString2rChar(axis$getNames()$toString()); return(er)})
            dim.list[[ensDimIndex]] <- list("Type" = axis$getAxisType()$toString(), "Units" = axis$getUnitsString(), "Values" = ens.values)
            names(dim.list)[ensDimIndex] <- "member" # axis$getDimensionsString()
      }
      if (grid$getTimeDimensionIndex() >= 0) {
            tDimIndex <- grid$getTimeDimensionIndex() + 1
            axis <- gcs$getTimeAxis()
            if (gcs$hasTimeAxis1D()) {
                  time.agg <- gcs$getTimeAxis1D()$getTimeResolution()$toString()
                  date.range <- axis$getCalendarDateRange()$toString()
                  # tdim.name <- axis$getDimensionsString()
            } else {
                  time.agg <- gcs$getTimeAxisForRun(0L)$getTimeResolution()$toString()
                  lastRunIndex <- as.integer(gcs$getRunTimeAxis()$getSize() - 1)
                  lastRun <- gcs$getTimeAxisForRun(lastRunIndex)
                  lastDateIndex <- as.integer(lastRun$getSize() - 1)    
                  time.end <- javaCalendarDate2rPOSIXlt(lastRun$getCalendarDate(lastDateIndex))
                  time.start <- javaCalendarDate2rPOSIXlt(gcs$getTimeAxisForRun(0L)$getCalendarDate(0L))
                  date.range <- range(time.start, time.end)
                  aux <- unlist(strsplit(axis$getDimensionsString(), split="\\s"))
                  #tdim.name <- aux[grep("time", aux)]
            }
            tdim.name <- "time"
            dim.list[[tDimIndex]] <- list("Type" = axis$getAxisType()$toString(), "TimeStep" = time.agg, "Units" = axis$getUnitsString(), "Date_range" = date.range)
            names(dim.list)[tDimIndex] <- tdim.name
      }
      if (grid$getZDimensionIndex() >= 0) {
            zDimIndex <- grid$getZDimensionIndex() + 1
            axis <- gcs$getVerticalAxis()
            gridLevels <- axis$getCoordValues()
            aux <- rep(NA, length(gridLevels))
            for (k in 1:length(gridLevels)) {
                  start <- rep(0L, length(gridShape))
                  count <- rep(1L, length(gridShape))
                  start[zDimIndex] <- as.integer(k - 1)
                  res <- tryCatch({grid$getVariable()$read(start, count)$getStorage()}, error = function(er) {err <- FALSE; return(err)})
                  aux[k] <- res
            }
            noLevelInd <- which(!aux)
            if (length(noLevelInd) > 0) {
                  gridLevels <- gridLevels[-noLevelInd]
            }
            dim.list[[zDimIndex]] <- list("Type" = axis$getAxisType()$toString(), "Units" = axis$getUnitsString(), "Values" = gridLevels)
            names(dim.list)[zDimIndex] <- "level" # axis$getDimensionsString()
      }
      if (grid$getXDimensionIndex() >= 0) {
            xDimIndex <- grid$getXDimensionIndex() + 1
            axis <- gcs$getXHorizAxis()
            lonAxisShape <- axis$getShape()        
            values <- axis$getCoordValues()
            if (length(lonAxisShape) == 2) {
                  values <- apply(t(matrix(values, ncol = lonAxisShape[1])), 2, min)
            }
            dim.list[[xDimIndex]] <- list("Type" = axis$getAxisType()$toString(), "Units" = axis$getUnitsString(), "Values" = values)
            names(dim.list)[xDimIndex] <- "lon" # axis$getDimensionsString()
      }
      if (grid$getYDimensionIndex() >= 0) {
            yDimIndex <- grid$getYDimensionIndex() + 1
            axis <- gcs$getYHorizAxis()
            latAxisShape <- axis$getShape()
            values <- axis$getCoordValues()
            if (length(latAxisShape) == 2) {
                  values <- apply(t(matrix(values, ncol = latAxisShape[1])), 1, min)
            }
            dim.list[[yDimIndex]] <- list("Type" = axis$getAxisType()$toString(), "Units" = axis$getUnitsString(), "Values" = values)
            names(dim.list)[yDimIndex] <- "lat" # axis$getDimensionsString()
      }
      return(dim.list)
}
# End
