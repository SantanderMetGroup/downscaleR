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
#' @param zRange A \sQuote{ucar.ma2.Range} or a null reference, as returned by
#'  \code{getVerticalLevelPars}
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
#' @references \url{https://www.unidata.ucar.edu/software/thredds/current/netcdf-java/v4.0/javadocAll/ucar/nc2/dt/grid/GeoGrid.html}
#' @author J. Bedia \email{joaquin.bedia@@gmail.com} and A. Cofin\~no.
#' @keywords internal
#' @export

makeSubset <- function(grid, timePars, zRange, latLon) {
      message("[", Sys.time(), "] Retrieving data subset ..." )
      gcs <- grid$getCoordinateSystem()
      dimNames <- rev(names(scanVarDimensions(grid)))
      aux.list <- list()
      for (i in 1:length(timePars$tRanges)) {
            dimNamesRef <- dimNames
            aux.list2 <- list()    
            for (j in 1:length(latLon$llbbox)) {
                  subSet <- grid$makeSubset(timePars$tRanges[[i]], zRange, latLon$llbbox[[j]], 1L, 1L, 1L)
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
            # Sub-routine for daily aggregation from 6h data
            if (!is.na(timePars$dailyAggr)) {
                  aux.list[[i]] <- toDD(do.call("abind", c(aux.list2, along = 1)), dimNamesRef, timePars$dailyAggr)
                  dimNamesRef <- attr(aux.list[[i]], "dimensions")
            } else {
                  aux.list[[i]] <- do.call("abind", c(aux.list2, along = 1))
            }
            aux.list2 <- NULL
      }
      mdArray <- do.call("abind", c(aux.list, along = grep("^time", dimNamesRef)))
      aux.list <- NULL
      if (any(dim(mdArray) == 1)) {
            dimNames <- dimNamesRef[-which(dim(mdArray) == 1)]    
            mdArray <- drop(mdArray)
      } else {
            dimNames <- dimNamesRef
      }
      mdArray <- unname(mdArray)
      attr(mdArray, "dimensions") <- dimNames
      return(mdArray)
}
# End