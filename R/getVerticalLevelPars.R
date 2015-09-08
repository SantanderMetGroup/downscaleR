#' Definition of vertical dimension slices
#'
#' Returns the selected level value (if any) and a suitable java structure. This is a subroutine
#' of \code{loadGridDataset}, whose output is passed to \code{makeSubset}.
#'
#' @param gcs An object of the java class \sQuote{GeoGrid})
#' @param level Vertical level. Passed by \code{loadGridDataset}, obtained via \code{findVerticalLevel}
#' 
#' @return A list with the level value and either a java Range or a java null reference
#' defining the index of the vertical axis (null if no vertical levels exist)
#' 
#' @details The function opens the grid and checks all possible level values.
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal
#' @import rJava


getVerticalLevelPars <- function(grid, level) {
      gcs <- grid$getCoordinateSystem()
      if (gcs$hasVerticalAxis()) {
            levels <- scanVarDimensions(grid)$level$Values
            if (is.null(level)) {
                  if (length(levels) == 1) {
                        level <- levels
                        if (gcs$getVerticalAxis()$findCoordElement(level) < 0) {
                              levelInd <- gcs$getVerticalAxis()$findCoordElement(0)
                        }
                  } else {
                        stop("Variable with vertical levels: '@level' following the variable name is required\nPossible values: ", paste(levels, collapse = ", "))
                  }
            } else {
                  if (level %in% levels) {
                        levelInd <- gcs$getVerticalAxis()$findCoordElement(level)
                  } else {
                        stop("Vertical level not found\nPossible values: ", paste(levels, collapse = ", "))
                  }
            }
            zRange <- .jnew("ucar/ma2/Range", levelInd, levelInd)
      } else {
            if (!is.null(level)) {
                  # warning("The variable selected is 2D: the '@level' specification was ignored")
                  level <- level
            }
            zRange <- .jnull()
      }
      return(list("level" = level, "zRange" = zRange))
}
# End
