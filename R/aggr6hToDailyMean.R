#' Performs the aggregation of 6 hourly data to daily mean
#' 
#' This is an auxiliary routine used in \code{\link{toDD}} within the \code{makeSubset} functions
#' 
#' @param x a vector of data (a time series actually)
#' @details The function assumes that we are loading 6-h data so that the length of x is always divisible by 4.
#' Otherwise it will raise an error. Because selections are based on seasons and years, there is no reason that this error
#'  ever occurs unless we are taking the previous forecast time for deaccumulation purposes (but so far there is not any
#'  de-accumulable 6h variable, and it is not likely that one is ever going to come across any...). 
#' @return The vector aggregated by the mean of each four elements
#' @references \url{http://stackoverflow.com/questions/18434305/apply-a-function-to-every-n-elements-of-a-vector}
#' @author J Bedia \email{joaquin.bedia@@gmail.com}
#' @export
#' @keywords internal

aggr6hToDailyMean <- function(x, aggr.fun) {
      stopifnot(length(x) %% 4 == 0L)
      if (aggr.fun == "sum") {
            dm <- tapply(x, rep(1 : (length(x) / 4), each = 4), sum)
      } else {
            dm <- tapply(x, rep(1 : (length(x) / 4), each = 4), mean)
      }
      return(dm)
}
# End


