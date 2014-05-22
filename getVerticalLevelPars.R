#' Definition of vertical dimension slices
#' 
#' Returns the selected level value (if any) and a suitable java structure. This is a subroutine
#' of \code{loadGridDataset}, whose output is passed to \code{makeSubset}. 
#'  
#' @param gcs An object of the java class \sQuote{GeoGrid})
#' @param level Vertical level. Passed by \code{loadGridDataset}
#' @return A list with the level value and either a java Range or a java null reference
#' defining the index of the vertical axis (null if no vertical levels exist)
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}

getVerticalLevelPars <- function(grid, level) {
    gcs <- grid$getCoordinateSystem()
    if (gcs$hasVerticalAxis()) {
        if (is.null(level)) {
            stop("Variable with vertical levels: Argument 'level' required.\nSee data inventory for valid argument values")
        }
        levelInd <- gcs$getVerticalAxis()$findCoordElement(level)
        if (levelInd < 0) {
            stop("Vertical level not found.\nCheck data inventory for valid vertical level values")
        }
        zRange <- .jnew("ucar/ma2/Range", levelInd, levelInd)
    } else {
        if (!is.null(level)) {
            warning("The variable selected has no vertical levels")
        }
        level <- NULL
        zRange <- .jnull()    
    }
    return(list("level" = level, "zRange" = zRange))
}
# End