#' Reads an arbitrary data slice
#' 
#' Reads an arbitrary data slice from a new GeoGrid that is a logical subset of
#' the original GeoGrid.
#' 
#' @param A grid of the java class \sQuote{ucar.nc2.dt.grid.GeoGrid}
#' @param tRanges A list of java ranges (\sQuote{ucar.ma2.Range} class) indicating
#'  the index positions in the time axis, as returned by \code{getTimeDomain}.
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
#' (either in X and/or in Y dimensions) must be adequately handled given taht the \sQuote{makeSubset} output 
#' has in this case a different shape than the final output after \dQuote{slicing}.
#' @references \url{https://www.unidata.ucar.edu/software/thredds/current/netcdf-java/v4.0/javadocAll/ucar/nc2/dt/grid/GeoGrid.html}
#' @author J. Bedia \email{joaquin.bedia@@gmail.com} and A. Cofin\~no.

# zRange <- levelPars$zRange

makeSubset <- function(grid, tRanges, zRange, latLon) {
    gcs <- grid$getCoordinateSystem()
    dimNames <- rev(names(scanVarDimensions(grid)))
    aux.list <- list()
    for (i in 1:length(tRanges)) {
        dimNamesRef <- dimNames
        aux.list2 <- list()    
        for (j in 1:length(latLon$llbbox)) {
            subSet <- grid$makeSubset(tRanges[[i]], zRange, latLon$llbbox[[j]], 1L, 1L, 1L)
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
#         aux.list[[i]] <- do.call("abind", c(aux.list2, along = grep(paste("^", gcs$getXHorizAxis()$getDimensionsString(), "$", sep = ""), dimNamesRef)[1]))
        aux.list[[i]] <- do.call("abind", c(aux.list2, along = 1))
        aux.list2 <- NULL
    }
    mdArray <- do.call("abind", c(aux.list, along = grep(gcs$getTimeAxis()$getDimensionsString(), dimNamesRef)))
    aux.list <- NULL
    if (any(dim(mdArray) == 1)) {
        dimNamesRef <- dimNamesRef[-which(dim(mdArray) == 1)]    
        mdArray <- unname(drop(mdArray))
    }
    attr(mdArray, "dimensions") <- dimNamesRef
    if (isTRUE(latLon$revLat)) {
        mdArray <- unname(revArrayLatDim(mdArray, grid))
        attr(mdArray, "dimensions") <- dimNamesRef
    }
    return(mdArray)
}
# End