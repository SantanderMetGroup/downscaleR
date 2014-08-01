#' @title Get regular grid definition 
#' @description Get the (regular) grid definition from an existing (gridded) dataset
#' @param gridData A grid data object coming from \code{\link{loadGridData}} or \code{\link{interpGridData}}
#'  or the function \code{\link[ecomsUDG.Raccess]{loadECOMS}} of package \pkg{ecomsUDG.Raccess}.
#' @author J Bedia \email{joaquin.bedia@@gmail.com}
#' @details The returned grid object inherits the attributes from the input \code{xyCoords} definition.
#' @export
#' @family loading.grid
#' 

getGrid <- function(gridData) {
      grid.x <- c(gridData$xyCoords$x[1], tail(gridData$xyCoords$x, 1))
      grid.y <- c(gridData$xyCoords$y[1], tail(gridData$xyCoords$y, 1))
      out <- list(x = grid.x, y = grid.y)
      attributes(out) <- attributes(gridData$xyCoords)
      attr(out,'resX')<-(tail(gridData$xyCoords$x, 1)-gridData$xyCoords$x[1])/(length(gridData$xyCoords$x)-1)
      attr(out,'resY')<-(tail(gridData$xyCoords$y, 1)-gridData$xyCoords$y[1])/(length(gridData$xyCoords$y)-1)
      return(out)
}
# End
