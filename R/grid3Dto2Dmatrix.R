#' @title Conversion of a 3D grid to a 2D matrix
#' @description Converts 3D grids of the form [time,lat,lon]
#'  to 2D matrices of the form [time, grid-point].
#' @param grid3D A 3-dimensional grid with longitude, latitude and time dimensions.
#' @return A 2-dimensional matrix with time in rows and grid-points in columns.
#' @details Global attributes of the output matrix are provided for metadata description and spatiotemporal collocation.
#' @author J. Bedia 
#' @export
# #' @seealso \code{\link{grid3Dto2Dmat}}, exported for the user.

grid3Dto2Dmatrix <- function(grid3D) {
      dimNames <- attr(grid3D$Data, "dimensions")
      if (is.null(dimNames)) stop("Undefined 'dimensions' attribute")
      if (!identical(c("time","lat","lon"), dimNames)) stop("Not a 3D grid")
      M <- array3Dto2Dmat(grid3D$Data)
      co <- if (!grepl("rotated", attr(grid3D$xyCoords, "projection"), ignore.case = TRUE)) {
            expand.grid(grid3D$xyCoords$y, grid3D$xyCoords$x)[2:1] 
      } else {
            expandIrregularCoords(grid3D$xyCoords)
      }
      attr(M, "varName") <- grid3D$Variable$varName
      attr(M, "level") <- grid3D$Variable$level
      attr(M, "units") <- attributes(grid3D$Variable)$"units"
      attr(M, "longname") <- attributes(grid3D$Variable)$"longname"
      attr(M, "Dates:start") <- grid3D$Dates$start     
      attr(M, "Dates:end") <-   grid3D$Dates$end   
      attr(M, "xyCoords:x") <- co[,1]      
      attr(M, "xyCoords:y") <- co[,2]
      attr(M, "xyCoords:projection") <- attr(grid3D$xyCoords, "projection")
      return(M)
}
# End


#' @title Expand irregular grid  coords
#' @description Expand the coordinates of irregular grids into a 2D matrix of x-y coordinates
#' @param irregularXYcoords And irregular xyCoords element of an irregular grid
#' @return a 2D matrix of x-y coordinates
#' @author J Bedia
#' @keywords internal

expandIrregularCoords <- function(irregularXYcoords) {
      lons <- irregularXYcoords$lon
      lats <- irregularXYcoords$lat
      coords <- expand.grid(irregularXYcoords$y, irregularXYcoords$x)[2:1]
      names(coords) <- c("x","y")
      co <- coords 
      for (i in 1:nrow(coords)) {
            row <- match(coords[i,2], irregularXYcoords$y)
            col <- match(coords[i,1], irregularXYcoords$x) 
            co[i,] <- c(irregularXYcoords$lon[row, col], irregularXYcoords$lat[row, col])
      }
      return(co)
}



#' #'@title Conversion 2D matrix into a 3D array
#' #'@description Converts a 2D matrix of the form [time, lonlat] to a 3D array of the form
#' #' [time,lat,lon], in this order. Mainly for PCA analysis and grid reconstruction.
#' #'@param mat2D A 2D matrix with time in rows and lonlat in columns, as returned 
#' #'by \code{\link{array3Dto2Dmat}} 
#' #'@param x unique X coordinates of the points, in ascending order
#' #'@param y As argument \code{x}, for the Y coordinates
#' #'@return A 3-dimensional array with the dimensions ordered: [time,lat,lon]
#' #'@importFrom abind abind
#' #'@details The function is the inverse of \code{\link{array3Dto2Dmat}} 
#' #'@author J. Bedia 
#' #'@keywords internal
#' #'@export
#' #'@seealso \code{\link{array3Dto2Dmat}}, which performs the inverse operation.
#' 
#' mat2Dto3Darray <- function(mat2D, x, y) {
#'       if (!all(x == sort(x)) | !all(y == sort(y))) {
#'             stop("Coordinates 'x' and 'y' must be given in ascending order")
#'       }
#'       aux.list <- lapply(1:nrow(mat2D), function(i) {
#'             t(matrix(mat2D[i, ], ncol = length(x), nrow = length(y)))
#'       })
#'       arr <- unname(do.call("abind", c(aux.list, along = -1)))
#'       aux.list <- NULL      
#'       arr <- aperm(arr, perm = c(1,3,2))
#'       attr(arr, "dimensions") <- c("time", "lat", "lon")
#'       return(arr)
#' }
#' # End
