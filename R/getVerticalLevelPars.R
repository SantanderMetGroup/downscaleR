#' Definition of vertical dimension slices
#'
#' Returns the selected level value (if any) and a suitable java structure. This is a subroutine
#' of \code{loadGridDataset}, whose output is passed to \code{makeSubset}.
#'
#' @param gcs An object of the java class \sQuote{GeoGrid})
#' @param level Vertical level. Passed by \code{loadGridDataset}, obtained via \code{findVerticalLevel}
#' @return A list with the level value and either a java Range or a java null reference
#' defining the index of the vertical axis (null if no vertical levels exist)
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal
#' @import rJava

getVerticalLevelPars <- function(grid, level) {
    gcs <- grid$getCoordinateSystem()
    if (gcs$hasVerticalAxis()) {
        if (is.null(level)) {
            stop("Variable with vertical levels: '@level' following the variable name is required")
        }
        levelInd <- gcs$getVerticalAxis()$findCoordElement(level)
        if (levelInd < 0) {
            stop("Vertical level not found")
        }
        zRange <- .jnew("ucar/ma2/Range", levelInd, levelInd)
    } else {
        if (!is.null(level)) {
            warning("The variable selected is 2D: the '@level' specification was ignored")
            level <- NULL
        }
        zRange <- .jnull()
    }
    return(list("level" = level, "zRange" = zRange))
}
# End
