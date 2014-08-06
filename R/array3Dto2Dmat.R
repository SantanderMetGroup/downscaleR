#'@title Coversion of a 3D array to a 2D matrix
#'@description Converts 3D arrays of the form [lon,lat,time] -not strictly in this order-,
#'  to 2D matrices of the form [time, grid-point], in this order. Mainly for PCA analysis.
#'@param array3D A 3-dimensional array with longitude, latitude and time dimensions
#'@return A 2-dimensional matrix with time in rows and grid-points in columns.
#'@details The function is intended to convert geo-grids to a convenient format for PCA-related analyses.
#' The columns are ordered in X and Y ascending order, with coordinate Y varying faster. Thus, column coordinates 
#' would be given by the expression: \emph{expand.grid(rev(gridData$xyCoords$y), gridData$xyCoords$x)[2:1]}. This
#' is the most convenient format in order to naturally fill a matrix with the adequate number of columns (longitudes) 
#' and rows (latitudes) given the vectorized value of the output at a given time (or after time-averaging via rowMeans).
#'@author J. Bedia \email{joaquin.bedia@@gmail.com}
#'@keywords internal
#'@export

array3Dto2Dmat <- function(array3D) {
      dimNames <- attr(array3D, "dimensions")
      lon.index <- grep("lon", dimNames)
      n.lon <- dim(array3D)[lon.index]
      lat.index <- grep("lat", dimNames)
      n.lat <- dim(array3D)[lat.index]
      aux.list <- list()
      indices <- rep(list(bquote()), length(dimNames))
      for (i in 1:n.lon) {
            indices[[lon.index]] <- i
            for (j in n.lat:1) { 
                  indices[[lat.index]] <- j
                  call <- as.call(c(list(as.name("["), quote(array3D)), indices))
                  aux.list[[length(aux.list) + 1]] <- eval(call)  
            }
      }
      M <- do.call("cbind", aux.list)
      return(M)
}
# End
