#' @title Get grid definition 
#' @description Get the grid definition from an existing (gridded) dataset
#' @param gridData A grid data object coming from \code{\link{loadGridData}} or \code{\link{interpGridData}}
#'  or the \code{ecomsUDG.Raccess} package function \code{\link[ecomsUDG.Raccess]{loadECOMS}}.
#' @author J Bedia \email{joaquin.bedia@@gmail.com}
#' @details The returned grid object inherits the attributes from the input \code{xyCoords} definition.
#' @export
#' @family loading.grid
#' 

getGrid <- function(gridData) {
      grid.x <- c(gridData$xyCoords$x[1], tail(gridData$xyCoords$x, 1), abs(gridData$xyCoords$x[2] - gridData$xyCoords$x[1]))
      grid.y <- c(gridData$xyCoords$y[1], tail(gridData$xyCoords$y, 1), abs(gridData$xyCoords$y[2] - gridData$xyCoords$y[1]))
      out <- list(grid.x, grid.y)
      attributes(out) <- attributes(gridData$xyCoords)
      return(out)
}
# End
