#' @title Plot a map of the mean value of a grid dataset
#' 
#' @description Plot the spatial mean of a gridded variable, or variables in the case of multi-fields.
#' 
#' @importFrom fields image.plot
#' @import maps 
#' @importFrom fields world
#' @param gridData A grid dataset
#' @param multi.member Should members be plotted sepparately (TRUE), or just a plot
#'  of the multi-member mean (FALSE, default)?. Ignored if the dataset has no members.
#' 
#' @return a plot of the mean/multi-member/multi-variable field with a world map superposed
#' @export
#' @details The function is a wrapper of the \code{\link[fields]{image.plot}} function
#' in package \pkg{fields}
#' @author J Bedia \email{joaquin.bedia@@gmail.com}
#' @note The function plots a simple temporal mean of the loaded object in the form of
#' a map. It does not handle other temporal aggregations. 
#' 
#' @examples
#' # A field
#' data(iberia_ncep_ta850)
#' plotMeanField(iberia_ncep_ta850)
#' # A multifield
#' data(iberia_ncep_hus850)
#' mf <- makeMultiField(iberia_ncep_ta850, iberia_ncep_hus850)
#' plotMeanField(mf)
#' # A multimember field
#' data(tasmax_forecast)
#' plotMeanField(tasmax_forecast) # multimember mean
#' plotMeanField(tasmax_forecast, multi.member = TRUE) # by members
#' # A multimember multifield
#' data(tasmin_forecast)
#' mm.mf <- makeMultiField(tasmax_forecast, tasmin_forecast)
#' plotMeanField(mm.mf) # Note: multi-member not supported in this case
#' 

plotMeanField <- function (gridData, multi.member = FALSE) {
      dimNames <- attr(gridData$Data, "dimensions")
      mar <- match(c("lon", "lat"), dimNames)
      if (length(mar) != 2) {
            stop("Not a rectangular spatial domain")
      }
      titles <- gridData$Variable$varName
      if (is.na(match("var", dimNames))) {
            if (("member" %in% dimNames) & isTRUE(multi.member)) {
                  titles <- gridData$Members
                  multiPlot(gridData, "member", titles)
            } else {
                  aux <- apply(gridData$Data, FUN = mean, MARGIN = mar, na.rm = TRUE)
                  image.plot(gridData$xyCoords$x, gridData$xyCoords$y, aux, xlab = "", ylab = "", asp = 1, horizontal = TRUE, cex.axis = .75)
                  title("")
                  mtext(titles)
                  world(add = TRUE)
            }
      } else {
            multiPlot(gridData, "var", titles)
      }
}
# End


       