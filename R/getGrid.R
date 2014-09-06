#' @title Get regular grid definition 
#' 
#' @description Get the (regular) grid definition from an existing (gridded) dataset
#' 
#' @param gridData A grid data object coming from \code{\link{loadGridData}} or \code{\link{interpGridData}}
#'  or the function \code{\link[ecomsUDG.Raccess]{loadECOMS}} of package \pkg{ecomsUDG.Raccess}.
#'  
#' @return A list of two named components, \code{x} and \code{y}, consisting of a vector of length two each one, defining
#' the x/y lower and upper bounds. The grid-cell resolution is given by the attributes \code{'resX'} and
#'  \code{'resY'} respectively.
#'  
#' 
#' @details In case of irregular grid definitions, the function forces the grid to regularity.
#' The returned grid object inherits the attributes from the input \code{xyCoords} definition.
#' 
#' @export
#' 
#' @family loading.grid
#' 
#' @author S. Herrera and J. Bedia \email{joaquin.bedia@@gmail.com}
#' 
#' @examples \donttest{
#' # Iberia domain
#' data(iberia_ncep_hus850)
#' getGrid(iberia_ncep_hus850)
#' # Europe domain
#' data(tasmin_forecast)
#' getGrid(tasmin_forecast)
#' plotMeanField(tasmin_forecast)
#' # Interpolate NCEP onto the System4 grid:
#' int <- interpGridData(iberia_ncep_hus850, getGrid(tasmin_forecast), "bilinear")
#' # Note the warnings because of the non-overlapping domain extents
#' plotMeanField(int)
#' # The other way round, using nearest neighbour interpolation:
#' int2 <- interpGridData(tasmin_forecast, getGrid(iberia_ncep_hus850))
#' plotMeanField(int2)
#' #' # In this case, the mismatch in domain extent occurs only in the longitudes (to the west)
#' }
#' 


getGrid <- function(gridData) {
      grid.x <- c(gridData$xyCoords$x[1], tail(gridData$xyCoords$x, 1))
      grid.y <- c(gridData$xyCoords$y[1], tail(gridData$xyCoords$y, 1))
      out <- list(x = grid.x, y = grid.y)
      attributes(out) <- attributes(gridData$xyCoords)
      if (!exists("resX", attributes(gridData$xyCoords))) {
            attr(out, "resX") <- (tail(gridData$xyCoords$x, 1) - gridData$xyCoords$x[1]) / (length(gridData$xyCoords$x) - 1)
      }
      if (!exists("resY", attributes(gridData$xyCoords))) {
            attr(out, "resY") <- (tail(gridData$xyCoords$y, 1) - gridData$xyCoords$y[1]) / (length(gridData$xyCoords$y) - 1)
      }
      return(out)
}
# End

