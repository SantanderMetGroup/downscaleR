#' Performs variable transformation
#' 
#' Performs variable transformation according to dictionary specifications.
#' Sub-routine of \code{loadGridDataset}.
#' 
#' @param dic Dictionary line for the variable, as returned by \code{dictionaryLookup}
#' @param timePars A list of time selection parameters, as returned by \code{getTimeDomain}
#' @param mdArray A n-dimensional array, as returned by \code{makeSubset}
#' @return a list with start/end dates and the transformed n-dimensional array of data.
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal

dictionaryTransform <- function(dic, grid, timePars, mdArray) {
      Data <- mdArray * dic$scale + dic$offset
	if (dic$deaccum != 0) {
	      dimName <- paste("^", grid$getCoordinateSystem()$getTimeAxis1D()$getName(), "$", sep = "")
            dimIndex <- grep(dimName, attr(mdArray, "dimensions"))
            t.ranges <- c(0, cumsum(sapply(1:length(timePars$tRanges), function(x) timePars$tRanges[[x]]$length())))
	      margin <- c(1:length(attr(mdArray, "dimensions")))[-dimIndex]
	      Data <- apply(Data, MARGIN = margin, FUN = deaccum, t.ranges)
      }
	return(Data)
}
# End
	 