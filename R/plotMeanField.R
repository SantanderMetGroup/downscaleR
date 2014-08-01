#' @title Plot a map of the mean value of a grid dataset
#' 
#' @description Plot the spatial mean of a gridded variable, or variables in the case of multi-fields.
#' 
#' @importFrom fields image.plot
#' @importFrom fields world
#' 
#' @param gridData A grid dataset
#' @param multi.member Should members be plotted sepparately (TRUE), or just a plot
#'  of the multi-member mean (FALSE, default)?. Ignored if the dataset has no members.
#' 
#' @return a plot of the mean/multi-member/multi-variable field with a world map superposed
#' @export
#' @details The function is a wrapper of the \code{\link[fields]{image.plot}} function
#' in package \pkg{fields}
#' @author J Bedia joaquin.bedia@@gmail.com
#' @note The function plots a simple temporal mean of the loaded object in the form of
#' a map. It does not handle other temporal aggregations. In case of multimember grid datasets,
#' It simply plots the multi-member mean.

plotMeanField <- function (gridData, multi.member = FALSE) {
      dimNames <- attr(gridData$Data, "dimensions")
      mar <- grep("lon|lat", dimNames)
      if (length(mar) != 2) {
            stop("Not a rectangular spatial domain")
      }
      if (is.na(match("var", dimNames))) {
            aux <- apply(gridData$Data, FUN = mean, MARGIN = mar)
            image.plot(gridData$xyCoords$x, gridData$xyCoords$y, aux, xlab = "", ylab = "", asp = 1, horizontal = TRUE)
            world(add = TRUE)
            if (("member" %in% dimNames) & isTRUE(multi.member)) {
                  multiPlot(gridData, "member")
            }
      } else {
            multiPlot(gridData, "var")
      }
}
# End
       