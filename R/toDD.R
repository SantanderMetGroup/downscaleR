#' Performs the 6h to 24h aggregation of variables
#' 
#' @param NDarray A N-dimensional array, as returned by \sQuote{sQuote}
#' @param dimNamesRefRef A vector of dimension names
#' @return A ND array aggregated by its time dimension
#' @details Because of the re-ordering of dimensions after using \code{apply}, the 
#' vector of dimension names is needed for re-arranging accordingly
#' @author J Bedia \email{joaquin.bedia@@gmail}

toDD <- function(NDarray, dimNamesRef, dailyAggr) {
      mar <- grep("^time", dimNamesRef, invert = TRUE)
      NDarray <- unname(apply(NDarray, MARGIN = mar, FUN = aggr6hToDailyMean, aggr.fun = dailyAggr))
      dimNamesRef <- c(grep("^time", dimNamesRef ,value = TRUE), dimNamesRef[-grep("^time", dimNamesRef)])
      attr(NDarray, "dimensions") <- dimNamesRef
      return(NDarray)
}
# End