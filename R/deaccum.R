#' Performs deaccumulation
#' 
#' Used for deaccumulation of precipitation in certain model data (e.g. System4). Sub-routine
#' of \code{dictionaryTransform}.
#' 
#' @param x a vector of (accumulated) data.
#' @param t.ranges a vector defining the start/end indices for each annual season, as returned by
#' \code{dictionaryTransform}. This is used to restart the deaccumulation at the beginning of each season.
#' @param dff Logical. Is deaccumulation performed from the first value of the time series?.
#'  Passed by dictionaryTransform,coming from getForecastTimeDomain.S4. See Details.
#' @return a vector of deaccumulated data
#' @details When leadMonth equals 0, there is not a previous day for starting deaccumulation, and therefore the
#'  first value of the time series is taken 'as is'. Otherwise, one value before the start has to be taken
#'  to preserve time series length (this is previously done by getForecastDomain.S4).
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal
#' @export

deaccum <- function(x, t.ranges, dff) {
      if (isTRUE(dff)) {
            l <- lapply(1:(length(t.ranges) - 1), function(i) {
                  c(x[t.ranges[i] + 1], diff(x[(t.ranges[i] + 1) : t.ranges[i + 1]]))
            })
      } else {
            l <- lapply(1 : (length(t.ranges) - 1), function(i) {
                  diff(x[(t.ranges[i] + 1) : t.ranges[i + 1]])
            })
      }
      return(unlist(l))
}
# End