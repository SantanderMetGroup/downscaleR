#' @title Plot a map of the mean value of a grid dataset
#' 
#' @description Plot the spatial mean of a gridded variable, or variables in the case of multi-fields.
#' 
#' @importFrom fields image.plot
#' @importFrom fields world
#' @importFrom graphics layout
#' @importFrom abind asub
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

plotMeanField <- function (gridData) {
      dimNames <- attr(gridData$Data, "dimensions")
      mar <- grep("lon|lat", dimNames)
      if (length(mar) != 2) {
            stop("Not a rectangular spatial domain")
      }
      if (is.na(match("var", dimNames))) {
            aux <- apply(gridData$Data, FUN = mean, MARGIN = mar)
            image.plot(gridData$xyCoords$x, gridData$xyCoords$y, aux, xlab = "", ylab = "", asp = 1)
            world(add = TRUE)
      } else {
            var.dim.index <- grep("var", dimNames, fixed = TRUE)
            n.vars <- length(gridData$Variable$varName)
            nrows <- ifelse(sqrt(n.vars) < round(sqrt(n.vars)), ceiling(sqrt(n.vars)), floor(sqrt(n.vars)))
            def.par <- par(no.readonly = TRUE)
            par(mfrow = dim(mat))
            # plot(2,2, ty = "n", axes = FALSE, xlab = "", ylab = "")
            for (i in 1:n.vars) {
                  aux <- asub(gridData$Data, idx = i, dims = var.dim.index)
                  mar <- grep("lon|lat", dimNames[-var.dim.index])
                  aux <- apply(aux, FUN = mean, MARGIN = mar)
                  # axes <- ifelse(i == 1, TRUE, FALSE)
                  image.plot(gridData$xyCoords$x, gridData$xyCoords$y, aux, xlab = "", ylab = "", asp = 1) #, axes = axes)
                  world(add = TRUE)
                  mtext(gridData$Variable$varName[i])
            }
            par(def.par)
            
      }
}
# End
            
            
       