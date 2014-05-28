#' Performs deaccumulation
#' 
#' Used for deaccumulation of precipitation in certain model data (e.g. System4). Sub-routine
#' of \code{dictionaryTransform}.
#' 
#' @param x a vector of (accumulated) data.
#' @param t.ranges a vector defining the start/end indices for each annual season, as returned ny
#' \code{dictionaryTransform}. This is used to restart the deaccumulation at the beginning of each season.
#' @return a vector of deaccumulated data
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}

deaccum <- function(x, t.ranges) {
    l <- lapply(1:(length(t.ranges) - 1), function(i) {
        c(x[t.ranges[i] + 1], diff(x[t.ranges[i] + 1 : (t.ranges[i + 1])]))
        })
    return(unlist(l))
}
# End