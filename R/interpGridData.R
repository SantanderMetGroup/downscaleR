#' @title Interpolate a dataset to a grid
#' 
#' @description Performs interpolation of a gridded dataset into a new user-defined grid using bilinear weights 
#' or nearest-neighbour methods.
#' 
#' @importFrom fields interp.surface.grid
#' 
#'  
#'  @param gridData A grid data object coming from \code{\link{loadGridData}} or the \pkg{ecomsUDG.Raccess} 
#'  package function \code{\link[ecomsUDG.Raccess]{loadECOMS}}.
#'  @param new.grid Definition of the new grid, in the form of a list with the x and y components, in thir order.
#'  Each component consists of a vector of length three with components \emph{from}, \emph{to} and \emph{by},
#'   in this order, similar as the arguments passed to the \code{\link[base]{seq}} function, giving the 
#'   westernmost, easternmost and grid cell width in the X axis and, in the same way,
#'   the southernmost, northernmost and grid cell resolution in the Y axis. See details.
#'  @param method Method for interpolation. Currently implemented methods are either \code{bilinear},
#'  for bilinear interpolation, and \code{nearest}, for nearest-neighbor interpolation.
#'  @return An interpolated object preserving the output structure of the input
#'   (See e.g. \code{\link{loadGridData}}) for details on the output structure. 
#'  @details  In case of default definition of either x, y or both grid coordinates, the default grid
#'  is calculated taking the corners of the current grid and assuming x and y resolutions equal to 
#'  the default \code{by} argument value in function \code{\link[base]{seq}}: \emph{by = ((to - from)/(length.out - 1))}.
#'  The bilinear interpolator is a wrapper of the \code{\link[fields]{interp.surface.grid}} function
#'  in package \pkg{fields}.
#'  The output has special attributes in the \code{xyCoords} element that indicate that the object
#'   has been interpolated. These attributes are \code{interpolation}, which indicates the method used and
#'   \code{resX} and \code{resY}, for the grid-cell resolutions in the X and Y axes respectively.
#'   It is also possible to pass the interpolator the grid of a previously existing grid dataset using the
#'   \code{\link{getGrid}} method.
#'  @note To avoid unnecessary NA values, the function will not extrapolate using a new grid outside the
#'  current extent of the dataset, returning an error message.
#'  @family loading.grid
#'  @author J. Bedia \email{joaquin.bedia@@gmail.com} and S. Herrera
#'  @export
#'  @examples \dontrun{
#' # This is the path to the package built-in NCEP dataset (assumes reading permission)
#'  ncep <- file.path(find.package("downscaleR"), "datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.ncml")
#' # Load air temperature at 1000 mb isobaric pressure level for boreal winter (DJF) 1991-2000
#' t1000.djf <- loadGridData(ncep, var = "ta@@100000", lonLim = c(-12,10), latLim = c(33,47), season = c(12,1,2), years = 1991:2000)
#' par(mfrow = c(2,1))
#' plotMeanField(t1000.djf)
#' # Bilinear interpolation to a smaller domain centered in Spain using a 0.5 degree resolution in bot X and Y axes
#' t1000.djf.05 <- interpGridData(t1000.djf, new.grid = list(x = c(-10,5,.5), y = c(36,44,.5)), method = "bilinear")
#' plotMeanField(t1000.djf.05)
#' par(mfrow=c(1,1))
#' # New attributes "interpolation", "resX" and "resY" indicate that the original data have been interpolated
#' str(t1000.djf.05$xyCoords)
#' }


interpGridData <- function(gridData, new.grid = list(x = NULL, y = NULL), method = c("bilinear", "nearest")) {
      x <- gridData$xyCoords$x
      y <- gridData$xyCoords$y
      # Definition of new grid
      if (is.null(new.grid$x)) {
        new.grid.x <- gridData$xyCoords$x
      }else{
        if (!exists("resX", where = attributes(new.grid))) {
          new.grid.x <- new.grid$x
        }else{
          if (length(new.grid$x) != 2 | new.grid$x[2] < new.grid$x[1]) {
            stop("Invalid grid definition in X")
          }
          if (new.grid$x[1] < floor(x[1]) | new.grid$x[2] > ceiling(tail(x, 1))) {
            warning("The new longitudes are outside the data extent")
          }
          new.grid.x <- do.call("seq", as.list(c(new.grid$x, attr(new.grid, 'resX'))))
        }
      }
      if (is.null(new.grid$y)) {
        new.grid.y <- gridData$xyCoords$y
      }else{
        if (!exists("resY", where = attributes(new.grid))) {
          new.grid.y <- new.grid$y
        }else{
          if (length(new.grid$y) != 2 | new.grid$y[2] < new.grid$y[1]) {
            stop("Invalid grid definition in Y")
          }
          if (new.grid$y[1] < floor(y[1]) | new.grid$y[2] > ceiling(tail(y, 1))) {
            warning("The new latitudes are outside the data extent")
          }
          new.grid.y <- do.call("seq", as.list(c(new.grid$y, attr(new.grid, 'resY'))))
        }
      }
#       if (exists("resX", where = attributes(new.grid))) {
#         if (is.null(new.grid$x)) {
#           new.grid.x <- seq(x[1], tail(x, 1))
#         } else {
#           if (length(new.grid$x) != 2 | new.grid$x[2] < new.grid$x[1]) {
#             stop("Invalid grid definition in X")
#           }
#           if (new.grid$x[1] < floor(x[1]) & new.grid$x[2] <= ceiling(tail(x, 1))) {
#             stop("The westernmost corner of the new grid is outside the data extent\n Minimum X accepted value: ", floor(x[1]))
#           }
#           if (new.grid$x[2] > ceiling(tail(x, 1)) & new.grid$x[1] >= floor(x[1])) {
#             stop("The easternmost corner of the new grid is outside the data extent\n Maximum X accepted value: ", ceiling(tail(x, 1)))
#           }
#           if (new.grid$x[2] > ceiling(tail(x, 1)) & new.grid$x[1] < floor(x[1])) {
#             stop("The new grid is outside the data extent\n Accepted X values in the range: [", floor(x[1]), ",", ceiling(tail(x, 1)), "]")
#           }
#           new.grid.x <- do.call("seq", as.list(c(new.grid$x, attr(new.grid, 'resX'))))
#         }
#       } else {
#         new.grid.x <- new.grid$x
#       }
#       if (exists("resY", where = attributes(new.grid))) {
#         if (is.null(new.grid$y)) {
#           new.grid.y <- seq(y[1], tail(y, 1))
#         } else {
#           if (length(new.grid$y) != 2 | new.grid$y[2] < new.grid$y[1]) {
#             stop("Invalid grid definition in Y")
#           }
#           if (new.grid$y[1] < floor(y[1]) & new.grid$y[2] <= ceiling(tail(y, 1))) {
#             stop("The southernmost corner of the new grid is outside the data extent\n Minimum Y accepted value: ", floor(y[1]))
#           }
#           if (new.grid$y[2] > ceiling(tail(y, 1)) & new.grid$y[1] >= floor(y[1])) {
#             stop("The northernmost corner of the new grid is outside the data extent\n Maximum Y accepted value: ", ceiling(tail(y, 1)))
#           }
#           if (new.grid$y[2] > ceiling(tail(y, 1)) & new.grid$y[1] < floor(y[1])) {
#             stop("The new grid is outside the data extent\n Accepted Y values in the range: [", floor(y[1]), ",", ceiling(tail(y, 1)), "]")
#           }
#           new.grid.y <- do.call("seq", as.list(c(new.grid$y, attr(new.grid, 'resY'))))
#         }
#       } else {
#         new.grid.y <- new.grid$y
#       }
      grid.list <- list("x" = new.grid.x, "y" = new.grid.y)
      new.grid.x <- NULL
      new.grid.y <- NULL
      if (any(grepl("member", attr(gridData$Data, "dimensions")))) {
            mem.ind <- grep("member", attr(gridData$Data, "dimensions"))
            n.members <- dim(gridData$Data)[mem.ind]
      } else {
            n.members <- NULL
      }
      time.ind <- grep("^time", attr(gridData$Data, "dimensions"))
      # Handles reverse ordering of lon-lat (i.e. lat-lon)
      ind.order <- match(c("lon", "lat"), attr(gridData$Data, "dimensions"))
      transpose <- ifelse(!identical(sort(ind.order), ind.order), TRUE, FALSE)
      if (is.null(n.members)) {
            message("[", Sys.time(), "] Performing ", method, " interpolation... may take a while")
            interp.list <- lapply(1:dim(gridData$Data)[time.ind], function(j) {
                  indices <- rep(list(bquote()), length(dim(gridData$Data)))
                  indices[[time.ind]] <- j
                  ind.call <- as.call(c(list(as.symbol("["), quote(gridData$Data)), indices))
                  z <- eval(ind.call)
                  if (isTRUE(transpose)) {
                        z <- t(z)      
                  }
                  if (method == "bilinear") {
                        int <- interp.surface.grid(list("x" = x, "y" = y, "z" = z), grid.list)$z
                  }
                  if (method == "nearest") {
                        int <- matrix(nrow = length(grid.list$x), ncol = length(grid.list$y))
                        for (k in 1:length(grid.list$x)) {
                              for (l in 1:length(grid.list$y)) {
                                    ind.x <- which.min(abs(x - grid.list$x[k]))
                                    ind.y <- which.min(abs(y - grid.list$y[l]))
                                    int[k,l] <- z[ind.x, ind.y]
                              }
                        }
                  }
                  z <- NULL
                  return(int)
            })
            gridData$Data <- unname(do.call("abind", c(interp.list, along = 3)))
            interp.list <- NULL
            attr(gridData$Data, "dimensions") <- c("lon", "lat", "time")
      } else {
            message("[", Sys.time(), "] Performing multi-member ", method, " interpolation... may take a while")
            aux.list <- list()
            for (i in 1:n.members) {
                  message("[", Sys.time(), "] Interpolating member ", i, " out of ", n.members)
                  interp.list <- lapply(1:dim(gridData$Data)[time.ind], function(j) {
                        indices <- rep(list(bquote()), length(dim(gridData$Data)))
                        indices[[time.ind]] <- j
                        indices[[mem.ind]] <- i
                        ind.call <- as.call(c(list(as.symbol("["), quote(gridData$Data)), indices))
                        z <- eval(ind.call)
                        if (isTRUE(transpose)) {
                              z <- t(z)      
                        }
                        if (method == "bilinear") {
                              int <- interp.surface.grid(list("x" = x, "y" = y, "z" = z), grid.list)$z
                        }
                        if (method == "nearest") {
                              int <- matrix(nrow = length(grid.list$x), ncol = length(grid.list$y))
                              for (k in 1:length(grid.list$x)) {
                                    for (l in 1:length(grid.list$y)) {
                                          ind.x <- which.min(abs(x - grid.list$x[k]))
                                          ind.y <- which.min(abs(y - grid.list$y[l]))
                                          int[k,l] <- z[ind.x, ind.y]
                                    }
                              }
                        }
                        z <- NULL
                        return(int)
                  })
                  aux.list[[i]] <- unname(do.call("abind", c(interp.list, along = 3)))
            }
            interp.list <- NULL
            gridData$Data <- unname(do.call("abind", c(aux.list, along = 4)))
            aux.list <- NULL
            attr(gridData$Data, "dimensions") <- c("lon", "lat", "time", "member")
      }      
      gridData$xyCoords[1:2] <- grid.list
      attr(gridData$xyCoords, "interpolation") <-  method
      attr(gridData$xyCoords, "resX") <- abs(grid.list$x[2] - grid.list$x[1])
      attr(gridData$xyCoords, "resY") <- abs(grid.list$y[2] - grid.list$y[1]) 
      message("[", Sys.time(), "] Done")
      return(gridData)
}
# End   
