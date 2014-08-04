#' @title Determine the geo-location parameters of an arbitrary user selection
#'
#' @description The function uses the \sQuote{GeoGrid} object and the parameters \code{lonLim}
#'  and \code{latLim} passed by \code{loadGridDataset} and calculates the corresponding
#'  index positions.
#'
#' @param grid Java class \sQuote{GeoGrid}
#' @param lonLim see \code{\link{loadGridDataset}}
#' @param latLim see \code{\link{loadGridDataset}}
#' 
#' @return A list with the following items
#' \itemize{
#'  \item{\code{llbbox}} A list of length 1 or two depending on whether the selected domain
#'  crosses or not the dateline and longitude units go from 0 to 360. See details.
#'  \item{\code{pointXYindex}} A vector of length two with the index positions of the
#'  selected XY coordinates -in this order- in case of point selections. See details.
#'  \item{\code{lonSlice}} The X coordinates of the domain selected
#'  \item{\code{latSlice}} The Y coordinates of the domain selected
#'  \item{\code{revLat}} Logical. Whether the order of latitudes should be reversed or
#'  not in order to map the data properly in the geographical space
#' }
#' 
#' @details In order to deal with the problem of dateline crossing, the selection is
#' partitioned into two, and the part of the domain with negative eastings is put in first
#' place for consistent spatial mapping.
#' The index position of lon and lat in the corresponding axes is returned
#' by \code{pointXYindex}, and is passed to the \sQuote{readDataSlice} method in
#'  \code{makeSubset}. For single point locations, this is a integer vector of length
#'   two defining these positions, while in the case of rectangular domains its value is
#'    c(-1L,-1L).
#'    
#' @note The function assumes that datasets have degrees-east and degress-north as units
#' of the corresponding X and Y axes.
#' 
#' @references \url{https://www.unidata.ucar.edu/software/thredds/current/netcdf-java/v4.0/javadocAll/ucar/nc2/dt/GridCoordSystem.html#getRangesFromLatLonRect\%28ucar.unidata.geoloc.LatLonRect\%29}
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com} and A. Cofin\~no
#' 
#' @keywords internal
#' 
#' @export
#' 
#' @import rJava

getLatLonDomain <- function(grid, lonLim, latLim) {
      if (any(lonLim > 180) | any(lonLim < -180) | any(latLim > 90) | any(latLim < -90)) {
            stop("Invalid geographical coordinates. Check 'lonLim' and/or 'latLim' argument values")
      }
      message("[", Sys.time(), "] Defining geo-location parameters")
      gcs <- grid$getCoordinateSystem()
      bboxDataset <- gcs$getLatLonBoundingBox()
      if (length(lonLim) == 1 | length(latLim) == 1) {
            pointXYpars <- findPointXYindex(lonLim, latLim, gcs)
            lonLim <- pointXYpars$lonLim
            latLim <- pointXYpars$latLim
            pointXYindex <- pointXYpars$pointXYindex
      } else {
            pointXYindex <- c(-1L, -1L)
      }
      if (is.null(lonLim)) {
            lonLim <- c(bboxDataset$getLonMin(), bboxDataset$getLonMax())
            if (any(lonLim > 180)) {
                  lonLim <- lonLim - 180
            }
      }
      if (is.null(latLim)) {
            latLim <- c(bboxDataset$getLatMin(), bboxDataset$getLatMax())
      }
      deltaLat <- latLim[2] - latLim[1]
      deltaLon <- lonLim[2] - lonLim[1]
      spec <- .jnew("java/lang/String", paste(latLim[1], lonLim[1], deltaLat, deltaLon, sep = ", "))
      bboxRequest <- .jnew("ucar/unidata/geoloc/LatLonRect", spec)
      llbbox <- list()
      if (bboxRequest$getLonMin() < 0 & bboxRequest$getLonMax() >= 0 & bboxDataset$crossDateline()) {
            spec1 <- .jnew("java/lang/String", paste(latLim[1], lonLim[1], deltaLat, 0 - lonLim[1], sep = ", "))
            spec2 <- .jnew("java/lang/String", paste(latLim[1], 0, deltaLat, lonLim[2], sep = ", "))
            llbbox[[1]] <- .jnew("ucar/unidata/geoloc/LatLonRect", spec1)
            llbbox[[2]] <- .jnew("ucar/unidata/geoloc/LatLonRect", spec2)
      } else {
            llbbox[[1]] <- .jnew("ucar/unidata/geoloc/LatLonRect", spec)
      }
      if (pointXYindex[1] >= 0) {
            aux <- grid$makeSubset(.jnull(), .jnull(), .jnull(), 1L, 1L, 1L)
            lonSlice <- aux$getCoordinateSystem()$getLonAxis()$getCoordValue(pointXYindex[1])
      } else {
            lonAux <- list()
            for (k in 1:length(llbbox)) {
                  aux <- grid$makeSubset(.jnull(), .jnull(), llbbox[[k]], 1L, 1L, 1L)
                  lonAxisShape <- aux$getCoordinateSystem()$getXHorizAxis()$getRank()
                  lonAux[[k]] <- aux$getCoordinateSystem()$getXHorizAxis()$getCoordValues()
                  if (length(lonAxisShape) > 1) {
                        lonAux[[k]] <- apply(t(matrix(lonAux[[k]], ncol = lonAxisShape[1])), 2, min)
                  }
            }
            lonSlice <- do.call("c", lonAux)
      }
      lonSlice[which(lonSlice > 180)] <- lonSlice[which(lonSlice > 180)] - 360
      lonSlice <- sort(lonSlice)
      aux <- grid$makeSubset(.jnull(), .jnull(), llbbox[[1]], 1L, 1L, 1L)
      revLat <- FALSE
      if (pointXYindex[2] >= 0) {
            latSlice <- aux$getCoordinateSystem()$getYHorizAxis()$getCoordValue(pointXYindex[2])
      } else {
            latSlice <- aux$getCoordinateSystem()$getYHorizAxis()$getCoordValues()
            latAxisShape <- aux$getCoordinateSystem()$getYHorizAxis()$getRank()
            if (length(latAxisShape) > 1) {
                  latSlice <- apply(t(matrix(latSlice, ncol = latAxisShape[1])), 1, min)
            }
            if (diff(latSlice)[1] < 0) {
                  latSlice <- rev(latSlice)
                  revLat <- TRUE
            }
      }
      return(list("llbbox" = llbbox, "pointXYindex" = pointXYindex, "xyCoords" = list("x" = lonSlice, "y" = latSlice), "revLat" = revLat))
}
# End
