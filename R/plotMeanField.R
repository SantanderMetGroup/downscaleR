#' Plot a map of the mean value of a grid dataset
#' 
#' @description Plot the mean of a loaded grid dataset. The funciton is intended for obtaining
#' a quick visual reference of the data loaded. In case of multi-member gridded datasets, it displays
#' the multi-member mean
#' 
#' @importFrom fields image.plot
#' @importFrom fields world
#' 
#' @param gridData A grid dataset
#' @return a plot of the mean field with a world map superposed
#' @export
#' @details The function is a wrapper of the \code{\link[fields]{image.plot}} function
#' in package \pkg{fields}
#' @author J Bedia joaquin.bedia@@gmail.com
#' @note The function plots a simple temporal mean of the loaded object in the form of
#' a map. It does not handle other temporal aggregations. In case of multimember grid datasets,
#' It simply plots the multi-member mean.
#' 
plotMeanField <- function (gridData) {
      dimNames <- attr(gridData$Data, "dimensions")
      mar <- grep("lon|lat", dimNames)
      if (length(mar) != 2) {
            stop("Not a rectangular spatial domain")
      }
      aux <- apply(gridData$Data, FUN = mean, MARGIN = mar)
      image.plot(gridData$xyCoords$x, gridData$xyCoords$y, aux, xlab = "", ylab = "", asp = 1)
      world(add = TRUE)
} 
# End
