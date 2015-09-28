#' @title Select an arbitrary subset from a field or multifield along one or more of its dimensions
#'
#' @description Create a new field/multifield that is a subset of the input field 
#' along the selected dimensions
#'
#' @param field The input field to be subset. This is either a field, as returned e.g. by \code{loadGridData}, a
#' multifield, as returned by \code{loadMultiField} or \code{makeMultiField}, or other types of multimember fields
#' (possibly multimember multifields) as returned e.g. by \code{\link[ecomsUDG.Raccess]{loadECOMS}}.
#' @param var Character vector indicating the variables(s) to be extracted. (Used for multifield subsetting). See details.
#' @param members An integer vector indicating \strong{the position} of the members to be subset.
#' @param years The years to be selected. Note that this can be either a continuous or discontinuous
#' series of years, the latter option often used in a cross-validation framework.
#'  See details for year-crossing seasons. Default to \code{NULL} (no subsetting is performed on the time dimension).
#' @param latLim Same as \code{lonLim} argument, but for latitude.
#' @param lonLim Vector of length = 2, with minimum and maximum longitude coordinates, in decimal degrees,
#'  of the bounding box defining the subset. For single-point subsets, a numeric value with the
#'  longitude coordinate. If \code{NULL} (default), no subsetting is performed on the longitude dimension
#' @return A new field object that is a logical subset of the input field along the specified dimensions.
#' @details
#' 
#' The attribute \code{subset} will be added to the different slots corresponding to the subset dimensions, taking
#' the value of the subroutine called in each case (e.g.: attribute subset will have the value \code{subsetSpatial}
#' in the xyCoords slot after spatial subsetting...).
#' 
#' \strong{Time slicing by years}
#' 
#' In case of year-crossing seasons (e.g. boreal winter (DJF), \code{season = c(12,1,2)}),
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
#'  \strong{Extracting fields from multifields}
#'  
#'  One or several variables from a multifield object can be extracted. Note that argument \code{var} is insensitive to the order of the variables, i.e.: variables
#' will be always returned in the same order they are in the original multifield.
#'  
#' @importFrom abind asub
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @export
#' @family subsetting
#' @examples
#' # Example 1 - Spatial / member subset
#' data(tasmax_forecast)
#' plotMeanField(tasmax_forecast, TRUE)
#' # Selection of a smaller domain over the Iberian Peninsula and members 3 and 7
#' sub <- subsetField(tasmax_forecast,
#'                    members = c(3,7),
#'                    lonLim = c(-10,5),
#'                    latLim = c(36,44))
#' plotMeanField(sub, multi.member = TRUE)
#' ## Example 2 - Subsetting a multimember multifield by variables
#' # Multimember multifield creation
#' data(tasmax_forecast)
#' data(tasmin_forecast)
#' data(tp_forecast)
#' mm.mf <- makeMultiField(tasmax_forecast, tasmin_forecast, tp_forecast)
#' plotMeanField(mm.mf)
#' # Extracting just minimum temperature
#' sub1 <- subsetField(mm.mf, var = "tasmin", members = 1:4)
#' plotMeanField(sub1, multi.member = TRUE)
#' # Extracting precipitation and maximum temperature
#' # (Note that the field variables are NOT re-ordered)
#' sub2 <- subsetField(mm.mf, var = c("tp", "tasmax"))
#' plotMeanField(sub2)


subsetField <- function(field, var = NULL, members = NULL, years = NULL, latLim = NULL, lonLim = NULL) {
      if (!is.null(var)) {
            field <- subsetVar(field, var)
      }
      if (!is.null(members)) {
            field <- subsetMembers(field, members)
      }
      if (!is.null(years)) {
            field <- subsetYears(field, years)
      }
      if (!is.null(lonLim) | !is.null(latLim)) {
            field <- subsetSpatial(field, lonLim, latLim)
      }
      return(field)
}
# End


#' Extract a field from a multifield object
#' 
#' Extracts a field from a multifield object. Multimember multifields are supported. Subroutine of subsetField.
#'
#' @param multiField Input multifield to be subset 
#' @param var Character vector indicating the variables(s) to be extracted
#' @return Either a (multimember)field or (multimember)multifield if one ore more variables
#' are selected respectively.
#' @details Argument \code{var} is insensitive to the order of the variables, i.e.: variables
#' will be always returned in the same order they are in the original multifield.
#' 
#' An attribute 'subset' with value 'subsetVar' is added to the Variable slot of the output subset.
#' 
#' @importFrom abind asub
#' @keywords internal
#' @export
#' @author J. Bedia (joaquin.bedia@@gmail.com)
#' @family subsetting

subsetVar <- function(multiField, var) {
      if (length(multiField$Variable$varName) == 1) {
            warning("Argument 'var' was ignored: Input field is not a multifield object")
            return(multiField)
      } 
      var.idx <- grep(paste0("^", var, "$", collapse = "|"), multiField$Variable$varName)
      if (length(var.idx) == 0) {
            stop("Variables indicated in argument 'var' not found")
      }
      if (length(var.idx) < length(var)) {
            stop("Some variables indicated in argument 'var' not found")
      }
      var.dim <- grep("var", attr(multiField$Data, "dimensions"))
      dimNames <- attr(multiField$Data, "dimensions")
      multiField$Data <- asub(multiField$Data, idx = var.idx, dims = var.dim, drop = TRUE)                  
      mf <- FALSE
      attr(multiField$Data, "dimensions") <- if (length(dim(multiField$Data)) == length(dimNames)) {
            mf <- TRUE
            dimNames
      } else {
            dimNames[-1]
      }
      multiField$Variable$varName <- multiField$Variable$varName[var.idx]
      multiField$Variable$level <- multiField$Variable$level[var.idx]
      attributes(multiField$Variable) <- lapply(attributes(multiField$Variable), "[", var.idx)
      multiField$Dates <- if (isTRUE(mf)) {
            multiField$Dates[var.idx]
      } else {
            multiField$Dates[[var.idx]]
      }
      attr(multiField$Variable, "subset") <- "subsetVar"
      return(multiField)
}
# End


#' Member subsets from a multimember field
#' 
#' Retrieves a field that is a logical subset of a multimember field along its 'member' dimension.
#'  Multimember multifields are supported. Subroutine of \code{\link{subsetField}}.
#'
#' @param mmField Input multimember field to be subset (possibly a multimember multifield).
#' @param members An integer vector indicating \strong{the position} of the members to be subset.
#' @return A field (or multifield) that is a logical subset of the input field along its 'member' dimension.
#' @details An attribute 'subset' with value 'subsetMembers' is added to the Members slot of the output subset.
#' @importFrom abind asub
#' @keywords internal
#' @export
#' @author J. Bedia (joaquin.bedia@@gmail.com)
#' @family subsetting

subsetMembers <- function(mmField, members = NULL) {
      dimNames <- attr(mmField$Data, "dimensions")
      if (length(grep("member", dimNames)) == 0) {
            warning("Argument 'members' was ignored: Input field is not a multimember field object")
            return(mmField)
      }      
      mem.dim <- grep("member", attr(mmField$Data, "dimensions"))
      if (!all(members %in% (1:dim(mmField$Data)[mem.dim]))) {
            stop("'members' dimension subscript out of bounds")
      }
      mmField$Data <- asub(mmField$Data, idx = members, dims = mem.dim, drop = TRUE)                  
      mf <- FALSE
      attr(mmField$Data, "dimensions") <- if (length(dim(mmField$Data)) == length(dimNames)) {
            mf <- TRUE
            dimNames
      } else {
            dimNames[-mem.dim]
      }
      mmField$Members <- mmField$Members[members]
      if (is.list(mmField$InitializationDates)) { # e.g. CFSv2 (members defined through lagged runtimes)
            mmField$InitializationDates <- mmField$InitializationDates[members]
      } 
      attr(mmField$Members, "subset") <- "subsetMembers"
      return(mmField)
}
# End


#' Year subsets from a multimember field
#' 
#' Retrieves a field that is a logical subset of a multimember field along its 'time' dimension,
#'  on a yearly basis. Multimember multifields are supported. Subroutine of \code{\link{subsetField}}.
#'
#' @param field Input field to be subset (possibly a multimember/multifield).
#' @param years An integer vector indicating the years to be subset.
#' @details An attribute 'subset' with value 'subsetYears' is added to the Dates slot of the output subset.
#' @return A field (or multifield) that is a logical subset of the input field along its 'time' dimension.
#' @importFrom abind asub
#' @keywords internal
#' @export
#' @author J. Bedia (joaquin.bedia@@gmail.com)
#' @family subsetting

subsetYears <- function(field, years = NULL) {
      dimNames <- attr(field$Data, "dimensions")
      all.years <- getYearsAsINDEX(field)
      aux.year.ind <- match(years, unique(all.years))
      if (length(intersect(years, all.years)) == 0) {
            stop("No valid years for subsetting. The argument \'years\' was ignored")
      }
      if (any(years < min(all.years) | years > max(all.years))) {
            stop("Some subset time boundaries outside the current field extent")
      }
      time.ind <- which(all.years %in% years)
      field$Data <- asub(field$Data, time.ind, grep("time", dimNames))
      attr(field$Data, "dimensions") <- dimNames
      # Verification Date adjustment
      field$Dates <- if (any(grepl("var", dimNames))) {
            lapply(1:length(field$Dates), function(i) {
                  lapply(field$Dates[[i]], function(x) x[time.ind])})
      } else {
            lapply(field$Dates, function(x) x[time.ind])
      }
      # Initialization time adjustment
      if ("member" %in% dimNames) {
            field$InitializationDates <- if (is.list(field$InitializationDates)) { # Lagged runtime config
                  lapply(field$InitializationDates, "[", aux.year.ind)      
            } else {
                  field$InitializationDates[aux.year.ind]
            }
      }
      attr(field$Dates, "subset") <- "subsetYears"
      return(field)
}
# End


#' Spatial subset from a field
#' 
#' Retrieves a field that is a logical subset of the input field along its 'lat' and 'lon' dimensions.
#'  Multimember multifields are supported. Subroutine of \code{\link{subsetField}}.
#'
#' @param field Input field to be subset (possibly a multimember multifield).
#' @param lonLim Vector of length = 2, with minimum and maximum longitude coordinates, in decimal degrees,
#'  of the bounding box defining the subset. For single-point subsets, a numeric value with the
#'  longitude coordinate. If \code{NULL} (default), no subsetting is performed on the longitude dimension
#' @param latLim Same as \code{lonLim} argument, but for latitude.
#' @details An attribute 'subset' with value 'subsetSpatial' is added to the xyCoords slot of the output subset.
#' @return A field (or multifield) that is a logical spatial subset of the input field.
#' @importFrom abind asub
#' @keywords internal
#' @export
#' @author J. Bedia (joaquin.bedia@@gmail.com)
#' @family subsetting
#' 
subsetSpatial <- function(field, lonLim = NULL, latLim = NULL) {
      dimNames <- attr(field$Data, "dimensions")
      if (!is.null(lonLim)) {
            if (!is.vector(lonLim) | length(lonLim) > 2) {
                  stop("Invalid longitudinal boundary definition")
            }
            lons <- getCoordinates(field)$x
            if (lonLim[1] < lons[1] | lonLim[1] > tail(lons, 1)) {
                  stop("Subset longitude boundaries outside the current field extent: \n(",
                       paste(getGrid(field)$x, collapse = ","), ")")
            }
            lon.ind <- which.min(abs(lons - lonLim[1]))
            if (length(lonLim) > 1) {
                  if (lonLim[2] < lons[1] | lonLim[2] > tail(lons, 1)) {
                        stop("Subset longitude boundaries outside the current field extent: \n(",
                             paste(getGrid(field)$x, collapse = ","), ")")
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
      }
      if (!is.null(latLim)) {
            if (!is.vector(latLim) | length(latLim) > 2) {
                  stop("Invalid latitudinal boundary definition")
            }
            lats <- getCoordinates(field)$y
            if (latLim[1] < lats[1] | latLim[1] > tail(lats, 1)) {
                  stop("Subset latitude boundaries outside the current field extent: \n(",
                       paste(getGrid(field)$y, collapse = ","), ")")
            }
            lat.ind <- which.min(abs(lats - latLim[1]))
            if (length(latLim) > 1) {
                  if (latLim[2] < lats[1] | latLim[2] > tail(lats, 1)) {
                        stop("Subset latitude boundaries outside the current field extent: \n(",
                             paste(getGrid(field)$y, collapse = ","), ")")
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
      }
      attr(field$xyCoords, "subset") <- "subsetSpatial"
      return(field)
}
# End


