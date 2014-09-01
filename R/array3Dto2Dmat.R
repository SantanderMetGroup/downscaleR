#'@title Conversion of a 3D array to a 2D matrix
#'@description Converts 3D arrays of the form [lon,lat,time] -not strictly in this order-,
#'  to 2D matrices of the form [time, grid-point], in this order. Mainly for PCA analysis.
#'@param array3D A 3-dimensional array with longitude, latitude and time dimensions
#'@return A 2-dimensional matrix with time in rows and grid-points in columns.
#'@details The function is intended to convert fields to a convenient format for PCA-related analyses.
#' The columns are ordered in X and Y ascending order, with coordinate Y varying faster. Thus, column coordinates 
#' would be given by the expression: \emph{expand.grid(gridData$xyCoords$y, gridData$xyCoords$x)[2:1]}. This
#' is the most convenient format in order to naturally fill a matrix with the adequate number of columns (longitudes) 
#' and rows (latitudes) given the vectorized value of the output at a given time (or after time-averaging via rowMeans).
#' The function is insensitive to the dimension ordering of the input data array 
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
            for (j in 1:n.lat) { 
                  indices[[lat.index]] <- j
                  call <- as.call(c(list(as.name("["), quote(array3D)), indices))
                  aux.list[[length(aux.list) + 1]] <- eval(call)  
            }
      }
      M <- do.call("cbind", aux.list)
      return(M)
}
# End


#'@title Conversion 2D matrix into a 3D array
#'
#'@description Converts a 2D matrix of the form [time, lonlat] to a 3D array of the form
#' [time,lat,lon], in this order. Mainly for PCA analysis and field reconstruction.
#'  
#'@param mat2D A 2D matrix with time in rows and lonlat in columns, as returned 
#'by \code{\link{array3Dto2Dmat}} 
#'@param x unique X coordinates of the points, in ascending order
#'@param y As argument \code{x}, for the Y coordinates
#'
#'@return A 3-dimensional array with the dimensions ordered: [time,lat,lon]
#'
#'@importFrom abind abind
#'
#'@details The function is the inverse of \code{\link{array3Dto2Dmat}} 
#'@author J. Bedia \email{joaquin.bedia@@gmail.com}
#'@keywords internal
#'@export
#'@seealso \code{\link{array3Dto2Dmat}}, which performs the inverse operation.

mat2Dto3Darray <- function(mat2D, x, y) {
      if (!all(x == sort(x)) | !all(y == sort(y))) {
            stop("Coordinates 'x' and 'y' must be given in ascending order")
      }
      aux.list <- lapply(1:nrow(mat2D), function(i) {
            t(matrix(mat2D[i, ], ncol = length(x), nrow = length(y)))
      })
      arr <- unname(do.call("abind", c(aux.list, along = -1)))
      aux.list <- NULL      
      arr <- aperm(arr, perm = c(1,3,2))
      attr(arr, "dimensions") <- c("time", "lat", "lon")
      return(arr)
}
# End

