##     subsetGrid.R Select an arbitrary subset from a grid
##
##     Copyright (C) 2016 Santander Meteorology Group (http://www.meteo.unican.es)
##
##     This program is free software: you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 3 of the License, or
##     (at your option) any later version.
## 
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
## 
##     You should have received a copy of the GNU General Public License
##     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Select an arbitrary subset from a grid or multigrid along one or more of its dimensions
#' @description Create a new grid/multigrid that is a subset of the input grid along the selected dimensions
#' @param grid The input grid to be subset. This is either a grid, as returned e.g. by \code{loadeR::loadGridData}, a
#' multigrid, as returned by \code{makeMultiGrid}, or other types of multimember grids
#' (possibly multimember grids) as returned e.g. by \code{loadeR.ECOMS::loadECOMS}.
#' @param var Character vector indicating the variables(s) to be extracted. (Used for multigrid subsetting). See details.
#' @param members An integer vector indicating \strong{the position} of the members to be subset.
#' @param runtime An integer vector indicating \strong{the position} of the runtimes to be subset.
#' @param years The years to be selected. Note that this can be either a continuous or discontinuous
#' series of years, the latter option often used in a cross-validation framework.
#'  See details for year-crossing seasons. Default to \code{NULL} (no subsetting is performed on the time dimension).
#' @param season An integer vector indicating the months to be subset. 
#' @param latLim Same as \code{lonLim} argument, but for latitude.
#' @param lonLim Vector of length = 2, with minimum and maximum longitude coordinates, in decimal degrees,
#'  of the bounding box defining the subset. For single-point subsets, a numeric value with the
#'  longitude coordinate. If \code{NULL} (default), no subsetting is performed on the longitude dimension
#' @param outside if TRUE subset coordinates outside the grid extent are allowed. Default is FALSE.
#' @param drop Logical (default is TRUE). Drop or keep dimensions of length 1.
#' @return A new grid object that is a logical subset of the input grid along the specified dimensions.
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
#'   for instance the \code{loadGridData} function, from package \pkg{loadeR}, with the exception that several checks are undertaken
#'   to ensure that the subset is actually within the current extent of the input grid. It is also possible to
#'   make single-point selections from a grid, just by specifying a single coordinate instead of a range
#'    as the argument value. For instance \code{lonLim = c(-10,10)} and \code{latLim = c(35,45)} indicate a
#'  rectangular window centered in the Iberian Peninsula), and single grid-cell values
#'  (for instance \code{lonLim = -3.21} and \code{latLim = 41.087} for retrieving the data in the closest grid
#'  point to the point coordinate -3.21E, 41.087N. In the last two cases, the function
#'  operates by finding the nearest (euclidean distance) grid-points to the coordinates introduced.
#'  
#'  \strong{Extracting grids from multigrids}
#'  
#'  One or several variables from a multigrid object can be extracted. Note that argument \code{var} is 
#'  insensitive to the order of the variables, i.e.: variables will be always returned in the same order
#'   they are in the original multigrid.
#'  
#' @importFrom abind asub
#' @author J. Bedia, M. Iturbide
#' @export
#' @family subsetting
#' @examples
#' # Example 1 - Spatial / member subset
#' data(tasmax_forecast)
#' plotMeanGrid(tasmax_forecast, TRUE)
#' # Selection of a smaller domain over the Iberian Peninsula and members 3 and 7
#' sub <- subsetGrid(tasmax_forecast,
#'                    members = c(3,7),
#'                    lonLim = c(-10,5),
#'                    latLim = c(36,44))
#' plotMeanGrid(sub, multi.member = TRUE)
#' ## Example 2 - Subsetting a multimember multigrid by variables
#' # Multimember multigrid creation
#' data(tasmax_forecast)
#' data(tasmin_forecast)
#' data(tp_forecast)
#' mm.mf <- makeMultiGrid(tasmax_forecast, tasmin_forecast, tp_forecast)
#' plotMeanGrid(mm.mf)
#' # Extracting just minimum temperature
#' sub1 <- subsetGrid(mm.mf, var = "tasmin", members = 1:4)
#' plotMeanGrid(sub1, multi.member = TRUE)
#' # Extracting precipitation and maximum temperature
#' # (Note that the grid variables are NOT re-ordered)
#' sub2 <- subsetGrid(mm.mf, var = c("tp", "tasmax"))
#' plotMeanGrid(sub2)


subsetGrid <- function(grid,
                       var = NULL,
                       runtime = NULL,
                       members = NULL,
                       years = NULL,
                       season = NULL,
                       lonLim = NULL,
                       latLim = NULL,
                       outside = FALSE,
                       drop = TRUE) {
      if (!is.null(var)) {
            grid <- subsetVar(grid, var, drop)
      }
      if (!is.null(runtime)) {
            grid <- subsetRuntime(grid, runtime, drop)
      }
      if (!is.null(members)) {
            grid <- subsetMembers(grid, members, drop)
      }
      if (!is.null(years)) {
            grid <- subsetYears(grid, years, drop)
      }
      if (!is.null(season)) {
            grid <- subsetSeason(grid, season, drop)
      }
      if (!is.null(lonLim) | !is.null(latLim)) {
            grid <- subsetSpatial(grid, lonLim, latLim, outside, drop)
      }
      return(grid)
}
# End


#' Extract a grid from a multigrid object
#' 
#' Extracts a grid from a multigrid object. Multimember multigrids are supported. Subroutine of subsetGrid.
#'
#' @param multiGrid Input multigrid to be subset 
#' @param var Character vector indicating the variables(s) to be extracted
#' @param drop Logical (default is TRUE). Drop or keep dimensions of length 1.
#' @return Either a (multimember)grid or (multimember)multigrid if one ore more variables
#' are selected respectively.
#' @details Argument \code{var} is insensitive to the order of the variables, i.e.: variables
#' will be always returned in the same order they are in the original multigrid.
#' 
#' An attribute 'subset' with value 'subsetVar' is added to the Variable slot of the output subset.
#' 
#' @importFrom abind asub
#' @keywords internal
#' @export
#' @author J. Bedia 
#' @family subsetting

subsetVar <- function(grid, var, drop) {
      if (length(grid$Variable$varName) == 1) {
            warning("Argument 'var' was ignored: Input grid is not a multigrid object")
            return(grid)
      } 
      var.idx <- grep(paste0("^", var, "$", collapse = "|"), grid$Variable$varName)
      if (length(var.idx) == 0) {
            stop("Variables indicated in argument 'var' not found")
      }
      if (length(var.idx) < length(var)) {
            stop("Some variables indicated in argument 'var' not found")
      }
      var.dim <- grep("var", attr(grid$Data, "dimensions"))
      dimNames <- attr(grid$Data, "dimensions")
      grid$Data <- asub(grid$Data, idx = var.idx, dims = var.dim, drop = drop)                  
      mf <- FALSE
      attr(grid$Data, "dimensions") <- if (length(dim(grid$Data)) == length(dimNames)) {
            mf <- TRUE
            dimNames
      } else {
            dimNames[-1]
      }
      grid$Variable$varName <- grid$Variable$varName[var.idx]
      grid$Variable$level <- grid$Variable$level[var.idx]
      attributes(grid$Variable)[-1] <- lapply(attributes(grid$Variable)[-1], "[", var.idx)
      grid$Dates <- if (isTRUE(mf)) {
            grid$Dates[var.idx]
      } else {
            grid$Dates[[var.idx]]
      }
      attr(grid$Variable, "subset") <- "subsetVar"
      return(grid)
}
# End


#' Member subsets from a multimember grid
#' 
#' Retrieves a grid that is a logical subset of a multimember grid along its 'member' dimension.
#'  Multimember multigrids are supported. Subroutine of \code{\link{subsetGrid}}.
#'
#' @param grid Input multimember grid to be subset (possibly a multimember multigrid).
#' @param members An integer vector indicating \strong{the position} of the members to be subset.
#' @param drop Logical (default is TRUE). Drop or keep dimensions of length 1.
#' @return A grid (or multigrid) that is a logical subset of the input grid along its 'member' dimension.
#' @details An attribute 'subset' with value 'subsetMembers' is added to the Members slot of the output subset.
#' @importFrom abind asub
#' @keywords internal
#' @export
#' @author J. Bedia 
#' @family subsetting

subsetMembers <- function(grid, members, drop) {
      dimNames <- attr(grid$Data, "dimensions")
      if (length(grep("member", dimNames)) == 0) {
            warning("Argument 'members' was ignored: Input grid is not a multimember grid object")
            return(grid)
      }      
      mem.dim <- grep("member", attr(grid$Data, "dimensions"))
      if (!all(members %in% (1:dim(grid$Data)[mem.dim]))) {
            stop("'members' dimension subscript out of bounds")
      }
      grid$Data <- asub(grid$Data, idx = members, dims = mem.dim, drop = drop)                  
      mf <- FALSE
      attr(grid$Data, "dimensions") <- if (length(dim(grid$Data)) == length(dimNames)) {
            mf <- TRUE
            dimNames
      } else {
            dimNames[-mem.dim]
      }
      grid$Members <- grid$Members[members]
      if (is.list(grid$InitializationDates)) { # e.g. CFSv2 (members defined through lagged runtimes)
            grid$InitializationDates <- grid$InitializationDates[members]
      } 
      if (!is.null(grid$Runtime)) attr(grid$Members, "subset") <- "subsetMembers"
      return(grid)
}
# End

#' Runtime subsets from a grid
#' 
#' Retrieves a grid that is a logical subset of a multiruntime grid along its 'runtime' dimension.
#'  Multiruntime multigrids are supported. Subroutine of \code{\link{subsetGrid}}.
#'
#' @param grid Input multiruntime grid to be subset.
#' @param runtime An integer vector indicating \strong{the position} of the runtimes to be subset.
#' @param drop Logical (default is TRUE). Drop or keep dimensions of length 1.
#' @return A grid (or multigrid) that is a logical subset of the input grid along its 'runtime' dimension.
#' @details An attribute 'subset' with value 'subsetRuntime' is added to the Runtime slot of the output subset.
#' @importFrom abind asub
#' @keywords internal
#' @export
#' @author M. Iturbide
#' @family subsetting

subsetRuntime <- function(grid, runtime, drop) {
      dimNames <- attr(grid$Data, "dimensions")
      if (length(grep("runtime", dimNames)) == 0) {
            warning("Argument 'runtime' was ignored: Input grid is not a multiruntime grid object")
            return(grid)
      }      
      run.dim <- grep("runtime", attr(grid$Data, "dimensions"))
      if (!all(runtime %in% (1:dim(grid$Data)[run.dim]))) {
            stop("'runtime' dimension subscript out of bounds")
      }
      grid$Data <- asub(grid$Data, idx = runtime, dims = run.dim, drop = drop)                  
      mf <- FALSE
      attr(grid$Data, "dimensions") <- if (length(dim(grid$Data)) == length(dimNames)) {
            mf <- TRUE
            dimNames
      } else {
            dimNames[-run.dim]
      }
      grid$Runtime <- grid$Members[runtime]
      if (is.list(grid$InitializationDates)) { # e.g. CFSv2 (members defined through lagged runtimes)
            grid$InitializationDates <- grid$InitializationDates[runtime]
      } 
      if (!is.null(grid$Runtime)) attr(grid$Runtime, "subset") <- "subsetRuntime"
      return(grid)
}
# End




#' Year subsets from a multimember grid
#' 
#' Retrieves a grid that is a logical subset of a multimember grid along its 'time' dimension,
#'  on a yearly basis. Multimember multigrids are supported. Subroutine of \code{\link{subsetGrid}}.
#'
#' @param grid Input grid to be subset (possibly a multimember/multigrid).
#' @param years An integer vector indicating the years to be subset.
#' @param drop Logical (default is TRUE). Drop or keep dimensions of length 1.
#' @details An attribute 'subset' with value 'subsetYears' is added to the Dates slot of the output subset.
#' @return A grid (or multigrid) that is a logical subset of the input grid along its 'time' dimension.
#' @importFrom abind asub
#' @keywords internal
#' @export
#' @author J. Bedia 
#' @family subsetting

subsetYears <- function(grid, years, drop) {
      dimNames <- attr(grid$Data, "dimensions")
      season <- getSeason(grid)
      all.years <- getYearsAsINDEX(grid)
      aux.year.ind <- match(years, unique(all.years))
      if (length(intersect(years, all.years)) == 0) {
            stop("No valid years for subsetting. The argument \'years\' was ignored")
      }
      if (any(years < min(all.years) | years > max(all.years))) {
            stop("Some subset time boundaries outside the current grid extent")
      }
      time.ind <- which(all.years %in% years)
      grid$Data <- asub(grid$Data, time.ind, grep("time", dimNames), drop = drop)
      attr(grid$Data, "dimensions") <- dimNames
      # Verification Date adjustment
      grid$Dates <- if (any(grepl("var", dimNames))) {
            lapply(1:length(grid$Dates), function(i) {
                  lapply(grid$Dates[[i]], function(x) x[time.ind])})
      } else {
            lapply(grid$Dates, function(x) x[time.ind])
      }
      # Initialization time adjustment
      if ("member" %in% dimNames) {
            grid$InitializationDates <- if (is.list(grid$InitializationDates)) { # Lagged runtime config
                  lapply(grid$InitializationDates, "[", aux.year.ind)      
            } else {
                  grid$InitializationDates[aux.year.ind]
            }
      }
      attr(grid$Dates, "subset") <- "subsetYears"
      attr(grid$Dates, "season") <- season
      return(grid)
}
# End


#' Spatial subset from a grid
#' 
#' Retrieves a grid that is a logical subset of the input grid along its 'lat' and 'lon' dimensions.
#'  Multimember multigrids are supported. Subroutine of \code{\link{subsetGrid}}.
#'
#' @param grid Input grid to be subset (possibly a multimember multigrid).
#' @param lonLim Vector of length = 2, with minimum and maximum longitude coordinates, in decimal degrees,
#'  of the bounding box defining the subset. For single-point subsets, a numeric value with the
#'  longitude coordinate. If \code{NULL} (default), no subsetting is performed on the longitude dimension
#' @param latLim Same as \code{lonLim} argument, but for latitude.
#' @param outside Logical. Default to \code{FALSE}. If \code{TRUE}, subset coordinates outside the grid extent 
#' are allowed. 
#' @param drop Logical (default to \code{TRUE}). Drop or keep dimensions of length 1.
#' @details An attribute \code{subset} with value \code{subsetSpatial} is added to the \code{xyCoords}
#' component of the output grid.
#' @return A grid (or multigrid) that is a logical spatial subset of the input grid.
#' @importFrom abind asub
#' @importFrom utils tail
#' @keywords internal
#' @export
#' @author J. Bedia 
#' @family subsetting
#' 
subsetSpatial <- function(grid, lonLim, latLim, outside, drop) {
      dimNames <- attr(grid$Data, "dimensions")
      if (!is.null(lonLim)) {
            if (!is.vector(lonLim) | length(lonLim) > 2) {
                  stop("Invalid longitudinal boundary definition")
            }
            lons <- getCoordinates(grid)$x
            if (lonLim[1] < lons[1] | lonLim[1] > tail(lons, 1)) {
                  if (outside == FALSE) {
                        stop("Subset longitude boundaries outside the current grid extent: \n(",
                        paste(getGrid(grid)$x, collapse = ","), ")")
                  } else {
                        warning("Subset longitude boundaries outside the current grid extent: \n(",
                          paste(getGrid(grid)$x, collapse = ","), ")")
                  }
            }
            lon.ind <- which.min(abs(lons - lonLim[1]))
            if (length(lonLim) > 1) {
                  if (lonLim[2] < lons[1] | lonLim[2] > tail(lons, 1)) {
                        if (outside == FALSE) {
                              stop("Subset longitude boundaries outside the current grid extent: \n(",
                              paste(getGrid(grid)$x, collapse = ","), ")")
                        } else {
                              warning("Subset longitude boundaries outside the current grid extent: \n(",
                              paste(getGrid(grid)$x, collapse = ","), ")")
                        }
                  }
                  lon2 <- which.min(abs(lons - lonLim[2]))
                  lon.ind <- lon.ind:lon2
                  grid$Data <- asub(grid$Data, lon.ind, grep("lon", dimNames))
                  attr(grid$Data, "dimensions") <- dimNames
            } else {
                  grid$Data <- asub(grid$Data, lon.ind, grep("lon", dimNames), drop = drop)
                  attr(grid$Data, "dimensions") <- dimNames[grep("lon", dimNames, invert = TRUE)]
                  dimNames <- attr(grid$Data, "dimensions")
            }
            grid$xyCoords$x <- grid$xyCoords$x[lon.ind]
      }
      if (!is.null(latLim)) {
            if (!is.vector(latLim) | length(latLim) > 2) {
                  stop("Invalid latitudinal boundary definition")
            }
            lats <- getCoordinates(grid)$y
            if (latLim[1] < lats[1] | latLim[1] > tail(lats, 1)) {
                  if (outside == FALSE) {
                        stop("Subset latitude boundaries outside the current grid extent: \n(",
                        paste(getGrid(grid)$y, collapse = ","), ")")
                  } else {
                        warning("Subset longitude boundaries outside the current grid extent: \n(",
                          paste(getGrid(grid)$x, collapse = ","), ")")
                  }
            }
            lat.ind <- which.min(abs(lats - latLim[1]))
            if (length(latLim) > 1) {
                  if (latLim[2] < lats[1] | latLim[2] > tail(lats, 1)) {
                        if (outside == FALSE) {
                              stop("Subset latitude boundaries outside the current grid extent: \n(",
                              paste(getGrid(grid)$y, collapse = ","), ")")
                        } else {
                              warning("Subset longitude boundaries outside the current grid extent: \n(",
                                paste(getGrid(grid)$x, collapse = ","), ")")
                        }
                  }
                  lat2 <- which.min(abs(lats - latLim[2]))
                  lat.ind <- lat.ind:lat2
                  grid$Data <- asub(grid$Data, lat.ind, grep("lat", dimNames))
                  attr(grid$Data, "dimensions") <- dimNames
            } else {
                  grid$Data <- asub(grid$Data, lat.ind, grep("lat", dimNames), drop = drop)
                  attr(grid$Data, "dimensions") <- dimNames[grep("lat", dimNames, invert = TRUE)]
                  dimNames <- attr(grid$Data, "dimensions")
            }
            grid$xyCoords$y <- grid$xyCoords$y[lat.ind]
      }
      attr(grid$xyCoords, "subset") <- "subsetSpatial"
      return(grid)
}
# End


#' Monthly subset of a grid
#' 
#' Retrieves a grid that is a logical subset of the input grid along its 'time' dimension,
#'  on a monthly basis. Multimember multigrids are supported. Subroutine of \code{\link{subsetGrid}}.
#'
#' @param grid Input grid to be subset (possibly a multimember/multigrid).
#' @param season An integer vector indicating the months to be subset.
#' @param drop Logical (default is TRUE). Drop or keep dimensions of length 1.
#' @details An attribute 'subset' with value 'subsetSeason' is added to the Dates slot of the output subset.
#' @return A grid (or multigrid) that is a logical subset of the input grid along its 'time' dimension.
#' @importFrom abind asub
#' @keywords internal
#' @export
#' @author J. Bedia 
#' @family subsetting

subsetSeason <- function(grid, season = NULL, drop = TRUE) {
      dimNames <- attr(grid$Data, "dimensions")
      season0 <- getSeason(grid)
      if (!all(season %in% season0)) stop("Month selection outside original season values")      
      mon <- if (any(grepl("var", dimNames))) {
            as.POSIXlt(grid$Dates[[1]]$start)$mon + 1
      } else {
            as.POSIXlt(grid$Dates$start)$mon + 1
      }
      time.ind <- which(mon %in% season)
      grid$Data <- asub(grid$Data, time.ind, grep("time", dimNames), drop = drop)
      attr(grid$Data, "dimensions") <- dimNames
      # Verification Date adjustment
      grid$Dates <- if (any(grepl("var", dimNames))) {
            lapply(1:length(grid$Dates), function(i) {
                  lapply(grid$Dates[[i]], function(x) x[time.ind])
            })
      } else {
            lapply(grid$Dates, function(x) x[time.ind])
      }
      attr(grid$Dates, "subset") <- "subsetSeason"
      return(grid)
}
# End


#' @title Select an arbitrary subset from a grid or multigrid along one of its dimensions
#' @description Create a new grid/multigrid that is a subset of the input grid 
#' along the selected dimension
#' @param grid The input grid to be subset. This is either a grid, as returned e.g. by \code{loadGridData} from package \pkg{loadeR} or a
#' multigrid, as returned by \code{makeMultiGrid}, or other types of multimember grids
#' (possibly multimember multigrids) as returned e.g. by \code{loadECOMS}, from package \pkg{loadeR.ECOMS}.
#' @param dimension Character vector indicating the dimension along which the positions indicated by the \code{indices} paraneter.
#' @param indices An integer vector indicating \strong{the positions} of the dimension to be extracted.
#' @param drop Logical (default is TRUE). Drop or keep dimensions of length 1.
#' @return A new grid object that is a logical subset of the input grid along the specified dimension.
#' @details
#' The attribute \code{subset} will be added taking the value of the \code{dimension} parameter.
#' @importFrom abind asub
#' @author J. Bedia and S. Herrera
#' @export
#' @family subsetting
#' @examples
#' # Example - Member subset
#' data(tasmax_forecast)
#' plotMeanGrid(tasmax_forecast, TRUE)
#' # Selection of a smaller domain over the Iberian Peninsula and members 3 and 7
#' sub <- subsetDimension(tasmax_forecast,
#'                    dimension = "member",
#'                    indices = c(1,3))
#' plotMeanGrid(sub, multi.member = TRUE)

subsetDimension <- function(grid, dimension = NULL, indices = NULL, drop = TRUE) {
      dimNames <- attr(grid$Data, "dimensions")
      if (!is.null(indices)) {
            grid$Data <- asub(grid$Data, indices, grep(dimension, dimNames), drop = drop)
            attr(grid$Data, "dimensions") <- dimNames
            if ("time" %in% dimension) {
                  grid$Dates$start <- grid$Dates$start[indices]
                  grid$Dates$end <- grid$Dates$end[indices]
            }
            if ("lon" %in% dimension) {
                  grid$xyCoords$x <- grid$xyCoords$x[indices]
            }
            if ("lat" %in% dimension) {
                  grid$xyCoords$y <- grid$xyCoords$y[indices]
            }
            if ("member" %in% dimension) {
                  grid$Members <- grid$Members[indices]
                  if (is.list(grid$InitializationDates)) { # e.g. CFSv2 (members defined through lagged runtimes)
                        grid$InitializationDates <- grid$InitializationDates[indices]
                  } 
            }
            attr(grid$Variable, "subset") <- dimension
      } else {
            warning("Argument 'indices' is NULL and no subsetting has been applied. The same input 'grid' is returned.")
      }
      return(grid)
}
# End


