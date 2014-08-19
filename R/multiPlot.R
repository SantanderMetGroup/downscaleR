#' @title Make multi-panel plots
#' 
#' @description Sub-routine of plotMeanField for dividing the graphical window into different subplots,
#' for multi-member or multi-variable displays
#' 
#' @param gridData a grid dataset as returned by any of the loading functions
#' @param name of the dimension used for splitting: either \code{var} or \code{member} for multi-predictor and
#' multi-member displays respectively
#' 
#' @return Prints the graphical display
#' 
#' @importFrom abind asub
#' @importFrom fields image.plot
#' @importFrom fields world
#' 
#' @keywords internal
#' @author J Bedia \email{joaquin.bedia@@gmail.com}


multiPlot <- function(gridData, split.dim.name, titles) {
      dimNames <- attr(gridData$Data, "dimensions")
      index <- grep(split.dim.name, dimNames, fixed = TRUE)
      n <- dim(gridData$Data)[index]
      nrows <- ifelse(sqrt(n) < round(sqrt(n)), ceiling(sqrt(n)), floor(sqrt(n)))
      mat <- matrix(1, ncol = ceiling(sqrt(n)), nrow = nrows)
      def.par <- par(no.readonly = TRUE)
      par(mfrow = dim(mat))
      for (i in 1:n) {
            aux <- asub(gridData$Data, idx = i, dims = index)
            mar <- match(c("lon", "lat"), dimNames[-index])
            aux <- apply(aux, mar, mean, na.rm = TRUE)
            image.plot(gridData$xyCoords$x, gridData$xyCoords$y, aux, xlab = "", ylab = "", asp = 1, horizontal = TRUE, cex.axis = .75) #, axes = axes)
            title("")
            mtext(titles[i])
            world(add = TRUE)
      }
      par(def.par)
}      
# End   
