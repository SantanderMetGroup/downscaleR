#' @title Reads an arbitrary data slice
#' 
#' @description Reads an arbitrary data slice from a new GeoGrid that is a logical subset of
#' the original GeoGrid.
#' 
#' @import rJava
#' @importFrom abind abind
#' 
#' @param A grid of the java class \sQuote{ucar.nc2.dt.grid.GeoGrid}
#' @param timePars A list of time parameters as returnde by \code{getTimeDomain}.
#' @param levelPars Vertical level parameters, as returned by \code{getVerticalLevelPars}.
#' The element \code{zRange} is a \sQuote{ucar.ma2.Range} or a null reference.
#' @param latLon A list of geospatial parameters, as returned by
#' \code{getLatLonDomain}.
#' @return A n-dimensional array with the selected subset data.
#' @note The process is somewhat tricky because R cannot adequately represent NDjavaArrays.
#' Several tests have shown that R (unlike MatLab) puts zeroes where they shouldn't be
#'  when using the \sQuote{copyToNDjavaArray} method. Hence, the data returned from
#'  \sQuote{readDataSlice} are converted to a 1D array (\sQuote{copyTo1DjavaArray} method),
#' and then put back in a n-dimensional R array. This is achieved by reversing
#' the GeoGrid dimension ordering. In addition, the particular case of single-coordinate selections
#' (either in X and/or in Y dimensions) must be adequately handled given that the \sQuote{makeSubset} output 
#' has in this case a different shape than the final output after \dQuote{slicing}.
#' @references \url{https://www.unidata.ucar.edu/software/thredds/v4.3/netcdf-java/v4.3/javadoc/ucar/nc2/dt/grid/GeoGrid.html}
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}, with contributions of A. Cofin\~no and S. Herrera.
#' @keywords internal
#' @export

makeSubset <- function(grid, timePars, levelPars, latLon) {
      message("[", Sys.time(), "] Retrieving data subset ..." )
      gcs <- grid$getCoordinateSystem()
      dimNames <- rev(names(scanVarDimensions(grid)))
      aux.list <- rep(list(bquote()), length(timePars$tRanges))
      do.aggr <- ifelse((timePars$aggr.d != "none") | (timePars$aggr.m != "none"), TRUE, FALSE)
      proj <- gcs$getProjection()
      for (i in 1:length(aux.list)) {
            dimNamesRef <- dimNames
            aux.list2 <- rep(list(bquote()), length(latLon$llbbox))
            for (j in 1:length(aux.list2)) {
                  if (!proj$isLatLon()){
                        subSet <- grid$makeSubset(levelPars$zRange, levelPars$zRange, timePars$tRanges[[i]], levelPars$zRange, latLon$latRanges, latLon$lonRanges)
                  }else{
                        subSet <- grid$makeSubset(timePars$tRanges[[i]], levelPars$zRange, latLon$llbbox[[j]], 1L, 1L, 1L)
                  }
                  shapeArray <- rev(subSet$getShape()) # Reversed!!
                  # shape of the output depending on spatial selection
                  if (latLon$pointXYindex[1] >= 0) {
                        rm.dim <- grep(gcs$getXHorizAxis()$getDimensionsString(), dimNamesRef, fixed = TRUE)
                        shapeArray <- shapeArray[-rm.dim]
                        dimNamesRef <- dimNamesRef[-rm.dim]
                  }
                  if (latLon$pointXYindex[2] >= 0) {
                        rm.dim <- grep(gcs$getYHorizAxis()$getDimensionsString(), dimNamesRef, fixed = TRUE)
                        shapeArray <- shapeArray[-rm.dim]
                        dimNamesRef <- dimNamesRef[-rm.dim]
                  }        
                  aux.list2[[j]] <- array(subSet$readDataSlice(-1L, -1L, latLon$pointXYindex[2], latLon$pointXYindex[1])$copyTo1DJavaArray(), dim = shapeArray)
            }
            aux.list[[i]] <- do.call("abind", c(aux.list2, along = 1))
            aux.list2 <- NULL
            # Daily aggregator 
            if (timePars$aggr.d != "none") {
                  aux.string <- paste(as.POSIXlt(timePars$dateSliceList[[i]])$mon, as.POSIXlt(timePars$dateSliceList[[i]])$mday, sep = "-")
                  aux.factor <- factor(aux.string, levels = unique(aux.string))
                  mar <- grep("^time", dimNamesRef, invert = TRUE)
                  aux.list[[i]] <- apply(aux.list[[i]], MARGIN = mar, FUN = function(x) {
                        tapply(x, INDEX = aux.factor, FUN = timePars$aggr.d, na.rm = TRUE)
                  })
                  dimNamesRef <- c("time", dimNamesRef[mar])
                  # Convert dates to daily:
                  nhours <- length(aux.factor) / nlevels(aux.factor)
                  timePars$dateSliceList[[i]] <- timePars$dateSliceList[[i]][seq(1, by = nhours, length.out = nlevels(aux.factor))]
            }
            # Monthly aggregator
            if (timePars$aggr.m != "none") {
                  mes <- as.POSIXlt(timePars$dateSliceList[[i]])$mon
                  day <- as.POSIXlt(timePars$dateSliceList[[i]])$mday
                  mar <- grep("^time", dimNamesRef, invert = TRUE)
                  aux.list[[i]] <- apply(aux.list[[i]], MARGIN = mar, FUN = function(x) {
                        tapply(x, INDEX = mes, FUN = timePars$aggr.m)
                  })
                  dimNamesRef <- if (length(unique(mes)) > 1) {
                        c("time", dimNamesRef[mar])
                  } else {
                        dimNamesRef[mar]
                  }
                  timePars$dateSliceList[[i]] <- timePars$dateSliceList[[i]][which(day == 1)]
            }
      }
      if (timePars$aggr.m != "none") {
            if (length(unique(mes)) > 1) {
                  mdArray <- do.call("abind", c(aux.list, along = grep("^time", dimNamesRef)))
            } else {
                  mdArray <- do.call("abind", c(aux.list, along = -1L))
                  dimNamesRef <- c("time", dimNamesRef)
            }
      } else {
            mdArray <- do.call("abind", c(aux.list, along = grep("^time", dimNamesRef)))
      }
      aux.list <- timePars$tRanges <- NULL
      if (any(dim(mdArray) == 1)) {
            dimNames <- dimNamesRef[-which(dim(mdArray) == 1)]    
            mdArray <- drop(mdArray)
      } else {
            dimNames <- dimNamesRef
      }
      mdArray <- unname(mdArray)
      attr(mdArray, "dimensions") <- dimNames
      # if (do.aggr) format <- ifelse(timePars$aggr.m == "none", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d")
      timePars$dateSliceList <- as.POSIXct(do.call("c", timePars$dateSliceList), tz = "GMT")#, format = format, tz = "GMT", usetz = TRUE)      
      return(list("timePars" = timePars, "mdArray" = mdArray))
}
# End