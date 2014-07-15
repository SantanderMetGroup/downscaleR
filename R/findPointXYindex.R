#' @title Define user-defined geolocation parameters
#' @description Finds the XY index position of single point selections
#' 
#' @param lonLim x domain definition
#' @param latLim y domain definition
#' @param gcs the grid coordinate system, formally a \sQuote{GridCoordSys} java object
#' @return a list with the new lonLim and latLim parameters (if changed) and the pointXYindex
#'  definition, consisting of a integer vector of length two.
#' @details The function applies the java method \sQuote{findXYindexFromCoord}, being the input
#' data the point coordinates in the native units. Thus, if longitudes are coded in the form [0-360] in the dataset,
#' the input lonLim coordinate [-180 - 180] is converted prior to passing it to the method.
#'  The resulting pointXYindex is ultimately passed to the java method \sQuote{readDataSlice} after subsetting.
#' The function also takes into account the possibility of selection of \dQuote{strips}, i.e.,
#' a selection across a whole or part of a meridian/parallel.
#' @author J Bedia \email{joaquin.bedia@@gmail.com}
#' @references \url{http://www.unidata.ucar.edu/software/thredds/v4.3/netcdf-java/v4.3/javadocAll/ucar/nc2/dt/grid/GridCoordSys.html#findXYindexFromLatLon\%28double,\%20double,\%20int[]\%29}
#' @keywords internal
#' @export
#' @import rJava


findPointXYindex <- function(lonLim, latLim, gcs)  {
      pointXYindex <- c(-1L, -1L)
      bboxDataset <- gcs$getLatLonBoundingBox()
      if (any(lonLim < 0) & bboxDataset$getLonMin() >= 0) {
            lon.aux <- sort(lonLim[which(lonLim < 0)] + 360)
      } else {
            if (!is.null(lonLim)) {
                  lon.aux <- lonLim
            } else {
                  lon.aux <- bboxDataset$getCenterLon()
            }
      }
      if (length(lonLim) == 1) {
            if (is.null(latLim)) {
                  lat.aux <- bboxDataset$getLatMin()
            } else {
                  lat.aux <- latLim[1]
            }
            pointXYindex[1] <- gcs$findXYindexFromCoord(lon.aux, lat.aux, .jnull())[1]
            if (pointXYindex[1] < 0) {
                  stop("Selected X point coordinate is out of range")
            }
            lonLim <- NULL   
      }
      if (length(latLim) == 1) {
            pointXYindex[2] <- gcs$findXYindexFromCoord(lon.aux[1], latLim, .jnull())[2]
            if (pointXYindex[2] < 0) {
                  stop("Selected Y point coordinate is out of range")
            }
            latLim <- NULL
      }
      return(list("lonLim" = lonLim, "latLim" = latLim, "pointXYindex" = pointXYindex))
}
# End