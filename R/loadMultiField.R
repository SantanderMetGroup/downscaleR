#' @title Load several gridded variables
#' 
#' @description Load a set of selected variables from the same datasets for a common 
#' spatio-temporal domain, typically to be used as predictors in a perfect-prog downscaling method.
#' 
#' @importFrom abind abind
#' 
#' @param
#' @param vars A character vector of length >= 2, specifying the names of the variables to be loaded.
#' @param 
#' @return
#' 
#' @details In essence, the function does as many calls to \code{\link{loadGridData}} as variables indicated in
#'  the \code{vars} argument, and then performs the interpolation and/or spatial checks to ensure the spatial consistency.
#' See \code{\link{loadGridData}} for details on input parameters for time, geolocation and homogenization parameters.
#' If any but not both of the components of the new.grid, the missing (NULL) one will be inherited from the grid of the
#'  first variable in \code{vars}
#' 
#' @export
#' 
#' @family loading
#' @family loading.grid
#' @family homogenization
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}


loadMultiField <- function(dataset, vars, dictionary = FALSE, lonLim = NULL, latLim = NULL, season = NULL, years = NULL, time = "none", new.grid = list(x = NULL, y = NULL), interp.method = "bilinear") {
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
                  x.aux <- apply(x.aux[ ,2:ncol(x.aux)], FUN = {function(x) {x - ref.x}}, MAR = 2)
                  y.aux <- apply(y.aux[ ,2:ncol(y.aux)], FUN = {function(x) {x - ref.y}}, MAR = 2)
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
                  message("[", Sys.time(), "] Regridding to 'new.grid' varible ", x, " (", vars[x], ") out of ", length(vars), " ...")
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
            lev <- var.list[[x]]$Variable$level
            if (is.null(lev)) {
                  lev <- NA
            }
            return(lev)
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
