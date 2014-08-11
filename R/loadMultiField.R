#' @title Load several gridded variables
#' 
#' @description Load a set of selected variables from the same datasets for a common 
#' spatio-temporal domain, typically to be used as predictors in a perfect-prog downscaling method.
#' 
#' @importFrom abind abind
#' 
#' @param dataset A character string indicating the database to be accessed. This is usually a path to a local file or a URL 
#' pointing to a netCDF or NcML file in the case of netCDF and/or gridded datasets. For station data in standard ASCII format,
#' this is the path to the directory the dataset lives in.
#' @param vars A character vector of length >= 2, specifying the names of the variables to be loaded. For variables with vertical
#'  levels, the vertical level is specified next to the variable name followed by the \dQuote{@@} symbol 
#'  (e.g. \code{var = "z@@700"} for geopotential heigth at 700 mb isobaric surface pressure level). It is also possible
#'   to enter the variable names as originally coded in the dataset to skip data homogenization (see nect argument).
#' @param dictionary Default to TRUE, indicating that a dictionary is used and the .dic file is stored in the same path than the
#' dataset. If the .dic file is stored elsewhere, then the argument is the full path to the .dic file (including the extension,
#' e.g.: \code{"/path/to/the/dictionary_file.dic"}). This is the case for instance when the dataset is stored in a remote URL,
#' and we have a locally stored dictionary for that particular dataset. If FALSE no variable homogenization takes place,
#' and the raw variable, as originally stored in the dataset, will be returned. See details for dictionary specification.
#' @param lonLim Vector of length = 2, with minimum and maximum longitude coordinates, in decimal degrees, of the bounding box selected.
#'  For single-point queries, a numeric value with the longitude coordinate. If \code{NULL} (default), the whole longitudinal range
#'   is selected (Note that this may lead to a large output object size).
#' @param latLim Same as \code{lonLim}, but for the selection of the latitudinal range.
#' @param season An integer vector specifying the desired season (in months, January = 1 ..., December = 12).
#'  Options include one to several (contiguous) months. Default to \code{NULL}, indicating a full year selection (same as \code{season = 1:12}).
#' @param years Optional vector of years to select. Default (\code{NULL}) to all available years. If the requested variable is static (e.g. orography)
#'  it will be ignored.  
#'  @param time A character vector indicating the temporal filtering/aggregation of the output data.
#'   Default to \code{"none"}, which returns the original time series as stored in the dataset. See \code{\link{loadGridData}}
#'   for further details.
#' \code{"00"}, \code{"06"}, \code{"12"} and \code{"18"}. If daily aggregated data are 
#' required use \code{"DD"}. If the requested variable is static (e.g. orography) it will be ignored. 
#' See details for time aggregation options
#' @param new.grid Definition of the new grid, in the form of a list with the x and y components, in this order.
#' If no \code{new.grid} is introduced, then the grid of the first variable in \code{vars} is taken. See details.
#' @param interp.method Interpolation method. Currently acceted values are \code{"bilinear"} and \code{"nearest"}. 
#' See \code{\link{interpGridData}} for details.
#'  
#' @return A list of components similar to the output of \code{\link{loadGridData}}, excepting for the fact that within the
#'  \code{Variable} slot, there is more than one variables in the \code{varName} and \code{level} elements. Also, the N-dimensional 
#'  array of the \code{Data} slot contains one additional dimension, corresponding to each of the variables loaded, following
#'   the order in the \code{Variable} slot info.
#' 
#' @details In essence, the function does as many calls to \code{\link{loadGridData}} as variables indicated in
#'  the \code{vars} argument, and then performs the interpolation and/or spatial checks to ensure the spatial consistency, via
#'  \code{\link{interpGridData}}. See \code{\link{loadGridData}} for details on input parameters for time, geolocation and homogenization parameters.
#' If only one (either x or y, but not both) of the components of the \code{new.grid} are provided, the missing one will be
#'  inherited from the grid of the first variable in \code{vars}, via \code{\link{getGrid}}.
#' 
#' @export
#' 
#' @family loading
#' @family loading.grid
#' @family homogenization
#' @family multifield
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}

loadMultiField <- function(dataset, vars, dictionary = TRUE, lonLim = NULL, latLim = NULL, season = NULL, years = NULL, time = "none", new.grid = list(x = NULL, y = NULL), interp.method = "bilinear") {
      if (length(vars) == 1) {
            stop("One single variable is not a multi-field.\nUse 'loadGridData' instead")
      }
      names(new.grid) <- c("x", "y")
      # loading vars
      var.list <- lapply (1:length(vars), function(x) {
            message("[", Sys.time(), "] Loading predictor ", x, " (", vars[x], ") out of ", length(vars))
            suppressMessages(loadGridData(dataset, vars[x], dictionary, lonLim, latLim, season, years, time))
      })
      # re-gridding
      if (is.null(new.grid$x) & !is.null(new.grid$y)) {
            new.grid$x <- getGrid(var.list[[1]])$x
      }
      if (is.null(new.grid$y) & !is.null(new.grid$x)) {
            new.grid$y <- getGrid(var.list[[1]])$y
      }
      if (is.null(new.grid$x) & is.null(new.grid$y)) {
            message("[", Sys.time(), "] Using the original grid as no 'new.grid' has been introduced")
            grid.list <- lapply(1:length(var.list), function(x) getGrid(var.list[[x]]))
            # Check for equality of grid definition
            aux <- do.call("c", grid.list)
            x.aux <- do.call("cbind", aux[seq(1, length(aux), 2)])
            y.aux <- do.call("cbind", aux[seq(2, length(aux), 2)])
            aux <- NULL
            ref.x <- x.aux[ ,1]
            ref.y <- y.aux[ ,1]
            if (ncol(x.aux) > 2) {
                  x.aux <- apply(x.aux[ ,2:ncol(x.aux)], 2, {function(x) {x - ref.x}})
                  y.aux <- apply(y.aux[ ,2:ncol(y.aux)], 2, {function(x) {x - ref.y}})
            } else {
                  x.aux <- x.aux[ ,2] - ref.x
                  y.aux <- y.aux[ ,2] - ref.y
            }
            ref.x <- NULL
            ref.y <- NULL
            if (sum(x.aux) != 0 | sum(y.aux) != 0) {
                  message("[", Sys.time(), "] Regridding to the grid of first variable (", vars[1], ") ...")
                  if (length(var.list) > 2) {
                        var.list <- lapply(2:length(var.list), function(x) {
                              var.list[[x]] <- suppressMessages(interpGridData(var.list[[x]], getGrid(var.list[[1]]), interp.method))
                        })
                  } else {
                        var.list[[2]] <- suppressMessages(interpGridData(var.list[[2]], getGrid(var.list[[1]]), interp.method))
                  }
            }
            x.aux <- NULL
            y.aux <- NULL
      } else {
            var.list <- lapply(1:length(var.list), function(x) {
                  message("[", Sys.time(), "] Regridding variable ", x, " (", vars[x], ") out of ", length(vars), " ...")
                  var.list[[x]] <- suppressMessages(interpGridData(var.list[[x]], new.grid, interp.method))
            })
      }
      # Variables
      varNames <- unlist(lapply(1:length(var.list), function(x) var.list[[x]]$Variable$varName))
      isStandard <- TRUE
      if (dictionary == FALSE) {
            isStandard <- FALSE
      } 
      level <- unlist(lapply(1:length(var.list), function(x) {
            return(ifelse(is.null(var.list[[x]]$Variable$level), NA, var.list[[x]]$Variable$level))
      }))
      Variable <- list("varName" = varNames, "isStandard" = isStandard, "level" = level)
      Dates <- var.list[[1]]$Dates
      xyCoords <- var.list[[1]]$xyCoords
      var.list <- lapply(1:length(var.list), function(x) {
            var.list[[x]] <- var.list[[x]]$Data
      })
      dimNames <- c("var", attr(var.list[[1]], "dimensions"))
      Data <- unname(do.call("abind", c(var.list, along = -1)))
      var.list <- NULL
      attr(Data, which = "dimensions") <- dimNames 
      message("[", Sys.time(), "] Done.")
      return(list("Variable" = Variable, "xyCoords" = xyCoords, "Data" = Data, "Dates" = Dates))
}
# End
