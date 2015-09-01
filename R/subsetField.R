#' @title Select an arbitrary spatial/temporal subset from a field or multifield
#'
#' @description Create a new field/multifield that is a subset of the input field 
#' along the space and/or time dimensions
#'
#' @param field The input data to be subset. This is either a field, as returned e.g. by \code{loadGridData}
#' or a Station dataset as returned by \code{loadStationData}.
#' @param lonLim Vector of length = 2, with minimum and maximum longitude coordinates, in decimal degrees,
#'  of the bounding box defining the subset. For single-point subsets, a numeric value with the
#'  longitude coordinate. If \code{NULL} (default), no subsetting is performed on the longitude dimension
#' @param latLim Same as \code{lonLim} argument, but for latitude.
#' @param years The years to be selected. Note that this can be either a continuous or discontinuous
#' series of years, the latter option often used in a cross-validation framework.
#'  See details for year-crossing seasons. Default to \code{NULL} (no subsetting is performed on the time dimension).
#' 
#' @return The same type of object as entered in the \code{field} argument representing a subset of its input
#' 
#' 
#' @details
#' 
#' \strong{Time slicing}
#' 
#' Time slicing is performed on a yearly basis (i.e., sub-yearly specifications are not allowed). In case of
#'  year-crossing seasons (e.g. boreal winter (DJF), \code{season = c(12,1,2)}),
#' the season is assigned to the years of January and February 
#' (i.e., winter of year 2000 corresponds to Dec 1999, Jan 2000 and Feb 2000). Thus, 
#' the \code{years} argument must be introduced accordingly (See e.g. \code{\link{getYearsAsINDEX}}
#' function for details).
#' 
#'  \strong{Spatial slicing}
#'  
#'  Spatial subset definition is done via the \code{lonLim} and \code{latLim} arguments, in the same way as
#'   for instance the \code{\link{loadGridData}} function, with the exception that several chechks are undertaken
#'   to ensure that the subset is actually within the current extent of the input field. It is also possible to
#'   make single-point selections from a field, just by specifying a single coordinate instead of a range
#'    as the argument value. For instance \code{lonLim=c(-10,10)} and \code{latLim=c(35,45)} indicates a
#'  rectangular window centered in the Iberian Peninsula), and single grid-cell values
#'  (for instance \code{lonLim=-3.21} and \code{latLim=41.087} for retrieving the data in the closest grid
#'  point to the point coordinate -3.21E, 41.087N. In the last two cases, the function
#'  operates by finding the nearest (euclidean distance) grid-points to the coordinates introduced.
#'  
#' @importFrom abind asub
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @export

subsetField <- function(field, lonLim = NULL, latLim = NULL, years = NULL) {
      dimNames <- attr(field$Data, "dimensions")
      if (!is.null(years)) {
            all.years <- getYearsAsINDEX(field)
            if (length(intersect(years, all.years)) == 0) {
                  stop("No valid years for subsetting. The argument \'years\' was ignored")
            }
            if (any(years < min(all.years) | years > max(all.years))) {
                  stop("Some subset time boundaries outside the current field extent")
            }
            time.ind <- which(all.years %in% years)
            field$Data <- asub(field$Data, time.ind, grep("time", dimNames))
            attr(field$Data, "dimensions") <- dimNames
            field$Dates <- if (any(grepl("var", dimNames))) {
                  lapply(1:length(field$Dates), function(i) {
                        lapply(field$Dates[[i]], function(x) x[time.ind])})
            } else {
                  lapply(field$Dates, function(x) x[time.ind])
            }
      }
      if (!is.null(lonLim)) {
            if (!is.vector(lonLim) | length(lonLim) > 2) {
                  stop("Invalid longitudinal boundary definition")
            }
            lons <- getCoordinates(field)$x
            if (lonLim[1] < lons[1] | lonLim[1] > tail(lons, 1)) {
                  stop("Subset longitude boundaries outside the current field extent: 
                       (", paste(getGrid(field)$x, collapse = ","), ")")
            }
            lon.ind <- which.min(abs(lons - lonLim[1]))
            if (length(lonLim) > 1) {
                  if (lonLim[2] < lons[1] | lonLim[2] > tail(lons, 1))
                        stop("Subset longitude boundaries outside the current field extent: 
                             (", paste(getGrid(field)$x, collapse = ","), ")")
            }
            lon2 <- which.min(abs(lons - lonLim[2]))
            lon.ind <- lon.ind:lon2
            field$Data <- asub(field$Data, lon.ind, grep("lon", dimNames))
            attr(field$Data, "dimensions") <- dimNames
      } else {
            field$Data <- asub(field$Data, lon.ind, grep("lon", dimNames), drop = TRUE)
            attr(field$Data, "dimensions") <- dimNames[grep("lon", dimNames, invert = TRUE)]
            dimNames <- attr(field$Data, "dimensions")
      }
      field$xyCoords$x <- field$xyCoords$x[lon.ind]
      if (!is.null(latLim)) {
            if (!is.vector(latLim) | length(latLim) > 2) {
                  stop("Invalid latitudinal boundary definition")
            }
            lats <- getCoordinates(field)$y
            if (latLim[1] < lats[1] | latLim[1] > tail(lats, 1)) {
                  stop("Subset latitude boundaries outside the current field extent: 
                 (", paste(getGrid(field)$x, collapse = ","), ")")
            }
            lat.ind <- which.min(abs(lats - latLim[1]))
            if (length(latLim) > 1) {
                  if (latLim[2] < lats[1] | latLim[2] > tail(lats, 1))
                        stop("Subset latitude boundaries outside the current field extent: 
                 (", paste(getGrid(field)$y, collapse = ","), ")")
            }
            lat2 <- which.min(abs(lats - latLim[2]))
            lat.ind <- lat.ind:lat2
            field$Data <- asub(field$Data, lat.ind, grep("lat", dimNames))
            attr(field$Data, "dimensions") <- dimNames
      } else {
            field$Data <- asub(field$Data, lat.ind, grep("lat", dimNames), drop = TRUE)
            attr(field$Data, "dimensions") <- dimNames[grep("lat", dimNames, invert = TRUE)]
            dimNames <- attr(field$Data, "dimensions")
      }
      field$xyCoords$y <- field$xyCoords$y[lat.ind]
      attr(field, "subset") <- "subsetField subset"
      return(field)
}
# End
