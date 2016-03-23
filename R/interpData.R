#' @title Interpolate a dataset
#' @description Interpolation of gridded datasets into a user-defined grid using nearest-neighbour or bilinear weights. 
#' @importFrom akima interp
#' @importFrom abind abind
#' @param obj A data object coming from \code{loadGridData}, \code{loadStationData} of package \pkg{loadeR}, or the \pkg{loadeR.ECOMS} 
#' package function \code{loadECOMS}.
#' @param new.Coordinates Definition of the new coordinates (grid or locations), in the form of a list with the x and y components, in thir order.
#' Each component consists of a vector of length three with components \emph{from}, \emph{to} and \emph{by},
#'  in this order, similar as the arguments passed to the \code{\link[base]{seq}} function, giving the 
#'  westernmost, easternmost and grid cell width in the X axis and, in the same way,
#'  the southernmost, northernmost and grid cell resolution in the Y axis. See details.
#' @param method Method for interpolation. Currently implemented methods are either \code{"bilinear"},
#' for bilinear interpolation, and \code{"nearest"}, for nearest-neighbor interpolation (default).
#' @return An interpolated object preserving the structure of the input
#' @details  In case of default definition of either x, y or both grid coordinates, the default grid
#' is calculated taking the corners of the current grid and assuming x and y resolutions equal to 
#' the default \code{by} argument value in function \code{\link[base]{seq}}: \emph{by = ((to - from)/(length.out - 1))}.
#' The bilinear interpolator uses the \code{\link[akima]{interp}} algorithm. 
#' The output has special attributes in the \code{xyCoords} element that indicate that the object
#'  has been interpolated. These attributes are \code{interpolation}, which indicates the method used and
#'  \code{resX} and \code{resY}, for the grid-cell resolutions in the X and Y axes respectively.
#'  It is also possible to pass the interpolator the grid of a previously existing grid dataset using the
#'  \code{\link{getGrid}} method.
#' @note To avoid unnecessary NA values, the function will not extrapolate using a new grid outside the
#' current extent of the dataset, returning an error message.
#' @author J. Bedia and S. Herrera
#' @export
#' @examples \dontrun{
#' # Load air temperature at 850 mb isobaric pressure level for boreal winter (DJF) 1991-2010
#' data(iberia_ncep_ta850)
#' par(mfrow = c(2,1))
#' plotMeanGrid(iberia_ncep_ta850)
#' # Bilinear interpolation to domain centered in Spain using a 0.5 degree resolution 
#' # in both X and Y axes
#' t <- interpData(iberia_ncep_ta850, new.Coordinates = list(x = c(-10,5,.5),
#'                                                           y = c(36,44,.5)),
#'                                    method = "bilinear")
#' plotMeanGrid(t)
#' par(mfrow=c(1,1))
#' # New attributes indicate that the data have been interpolated:
#' attributes(t$xyCoords)
#' }

interpData <- function(obj, new.Coordinates = list(x = NULL, y = NULL), method = c("nearest", "bilinear"), parallel, max.ncores, ncores) {
      method <- match.arg(method, choices = c("nearest", "bilinear"))
      parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
      if (any(attr(obj$Data, "dimensions") == "station")) {
            x <- as.numeric(obj$xyCoords[,1])
            y <- as.numeric(obj$xyCoords[,2])
      } else {
            x <- obj$xyCoords$x
            y <- obj$xyCoords$y
            if (!is.null(attr(obj$xyCoords, which = "projection"))) {
                  if (attr(obj$xyCoords, which = "projection") == "RotatedPole") {
                        x <- obj$xyCoords$lon
                        y <- obj$xyCoords$lat
                  }
            }
            if (is.vector(x) & is.vector(y)) {
                  x <- t(list(x = outer(obj$xyCoords$y*0,obj$xyCoords$x,FUN = "+"),
                              y = outer(obj$xyCoords$y,obj$xyCoords$x*0,FUN = "+"))$x)
                  y <- t(list(x = outer(obj$xyCoords$y*0,obj$xyCoords$x,FUN = "+"),
                              y = outer(obj$xyCoords$y,obj$xyCoords$x*0,FUN = "+"))$y)
            }
      }
      if (is.null(new.Coordinates)) {
            new.Coordinates <- getGrid(obj)
      } else {
            if (is.null(new.Coordinates$x)) {
                  new.Coordinates$x <- x
            } else if (exists("resX", where = attributes(new.Coordinates))) {
                  if (length(new.Coordinates$x) != 2 | new.Coordinates$x[2] < new.Coordinates$x[1]) {
                        stop("Invalid grid definition in X")
                  }
                  if ((max(c(new.Coordinates$x[1],new.Coordinates$x[2])) < min(x)) | (min(c(new.Coordinates$x[1],new.Coordinates$x[2])) > max(x))) {
                        stop("The input and output grids do not overlap\nCheck the input and output grid definitions")
                  }
                  if (new.Coordinates$x[1] < floor(min(x)) | new.Coordinates$x[2] > ceiling(max(x))) {
                        warning("The new longitudes are outside the data extent")
                  }
                  new.Coordinates$x <- do.call("seq", as.list(c(new.Coordinates$x, attr(new.Coordinates, 'resX'))))
            } else if (length(new.Coordinates$x) == 3) {
                  if (new.Coordinates$x[2] < new.Coordinates$x[1]) {
                        stop("Invalid grid definition in X")
                  }
                  if ((max(c(new.Coordinates$x[1],new.Coordinates$x[2])) < min(x)) | (min(c(new.Coordinates$x[1],new.Coordinates$x[2])) > max(x))) {
                        stop("The input and output grids do not overlap\nCheck the input and output grid definitions")
                  }
                  if (new.Coordinates$x[1] < floor(min(x)) | new.Coordinates$x[2] > ceiling(max(x))) {
                        warning("The new longitudes are outside the data extent")
                  }
                  if ((new.Coordinates$x[2] > new.Coordinates$x[1]) & (abs(new.Coordinates$x[3]) < abs(new.Coordinates$x[2] - new.Coordinates$x[1]))) {
                        new.Coordinates$x <- seq(from = new.Coordinates$x[1], to = new.Coordinates$x[2], by = new.Coordinates$x[3])
                  }
            }
            if (is.null(new.Coordinates$y)) {
                  new.Coordinates$y <- y
            } else if (exists("resY", where = attributes(new.Coordinates))) {
                  if (length(new.Coordinates$y) != 2 | new.Coordinates$y[2] < new.Coordinates$y[1]) {
                        stop("Invalid grid definition in Y")
                  }
                  if ((max(c(new.Coordinates$y[1],new.Coordinates$y[2])) < min(y)) | (min(c(new.Coordinates$y[1],new.Coordinates$y[2])) > max(y))) {
                        stop("The input and output grids do not overlap\nCheck the input and output grid definitions")
                  }
                  if (new.Coordinates$y[1] < floor(min(y)) | new.Coordinates$y[2] > ceiling(max(y))) {
                        warning("The new latitudes are outside the data extent")
                  }
                  new.Coordinates$y <- do.call("seq", as.list(c(new.Coordinates$y, attr(new.Coordinates, 'resY'))))
            } else if (length(new.Coordinates$y) == 3) {
                  if (new.Coordinates$y[2] < new.Coordinates$y[1]) {
                        stop("Invalid grid definition in Y")
                  }
                  if ((max(c(new.Coordinates$y[1],new.Coordinates$y[2])) < min(y)) | (min(c(new.Coordinates$y[1],new.Coordinates$y[2])) > max(y))) {
                        stop("The input and output grids do not overlap\nCheck the input and output grid definitions")
                  }
                  if (new.Coordinates$y[1] < floor(min(y)) | new.Coordinates$y[2] > ceiling(max(y))) {
                        warning("The new latitudes are outside the data extent")
                  }
                  if ((new.Coordinates$y[2] > new.Coordinates$y[1]) & (abs(new.Coordinates$y[3]) < abs(new.Coordinates$y[2] - new.Coordinates$y[1]))) {
                        new.Coordinates$y <- seq(from = new.Coordinates$y[1], to = new.Coordinates$y[2], by = new.Coordinates$y[3])
                  }
            }
      }
      if (any(grepl("member", attr(obj$Data, "dimensions")))) {
            mem.ind <- grep("member", attr(obj$Data, "dimensions"))
            n.members <- dim(obj$Data)[mem.ind]
      } else {
            n.members <- NULL
      }
      time.ind <- grep("^time", attr(obj$Data, "dimensions"))
      if (!any(attr(obj$Data, "dimensions") == "station")) {
            # Handles reverse ordering of lon-lat (i.e. lat-lon)
            ind.order <- match(c("lon", "lat"), attr(obj$Data, "dimensions"))
            transpose <- ifelse(!identical(sort(ind.order), ind.order), TRUE, FALSE)
      } else {
            transpose <- FALSE
      }
      if (is.null(n.members)) {
            message("[", Sys.time(), "] Performing ", method, " interpolation... may take a while")
            # To reduce the copy & paste code we use a kind of wrapper
            # function for lapply 
            if (parallel.pars$hasparallel) {
                apply_fun    <- function(...) {
                    parallel::parLapply(cl = parallel.pars$cl, ...)
                }  
                on.exit(parallel::stopCluster(parallel.pars$cl))
            } else {
                apply_fun <- lapply
            }
            interp.list <- apply_fun(1:dim(obj$Data)[time.ind], function(j) {
                  indices <- rep(list(bquote()), length(dim(obj$Data)))
                  indices[[time.ind]] <- j
                  ind.call <- as.call(c(list(as.symbol("["), quote(obj$Data)), indices))
                  z <- eval(ind.call)
                  if (isTRUE(transpose)) {
                        z <- t(z)
                  }
                  indNoNA <- which(!is.nan(z) | !is.na(z))
                  if (any(attr(new.Coordinates,"type") == "location")) {
                        int <- array(data = NA, dim = length(new.Coordinates$x))
                  }else{
                        int <- matrix(data = NA, nrow = length(new.Coordinates$x), ncol = length(new.Coordinates$y))
                  }
                  if (method == "bilinear" & length(indNoNA) != 0) {
                        int <- interp(x = x[indNoNA], y = y[indNoNA], z[indNoNA], xo = new.Coordinates$x, yo = new.Coordinates$y,
                                      linear = TRUE, extrap = FALSE, duplicate = "error", dupfun = NULL, ncp = NULL,
                                      nx = length(new.Coordinates$x), ny = length(new.Coordinates$y))
                        if (any(attr(new.Coordinates,"type") == "location")) {
                              int <- int$z[c(0:(length(new.Coordinates$y) - 1)) * length(new.Coordinates$y) + c(1:length(new.Coordinates$y))]
                        } else {
                              int <- int$z
                        }
                  }
                  if (method == "nearest" & length(indNoNA) != 0) {
                        if (any(attr(new.Coordinates,"type") == "location")) {
                              for (k in 1:length(new.Coordinates$x)) {
                                    if (!any(attr(obj$Data, "dimensions") == "station")) {
                                          ind.x <- arrayInd(which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[k]) ^ 2)), dim(x))[1]
                                          ind.y <- arrayInd(which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[k]) ^ 2)), dim(x))[2]
                                          int[k] <- z[ind.x, ind.y]
                                    } else {
                                          ind.x <- which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[k]) ^ 2))
                                          int[k] <- z[ind.x]
                                    }
                              }
                        } else {
                              int <- matrix(nrow = length(new.Coordinates$x), ncol = length(new.Coordinates$y))
                              for (k in 1:length(new.Coordinates$x)) {
                                    for (l in 1:length(new.Coordinates$y)) {
                                          if (!any(attr(obj$Data, "dimensions") == "station")) {
                                                ind.x <- arrayInd(which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[l]) ^ 2)), dim(x))[1]
                                                ind.y <- arrayInd(which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[l]) ^ 2)), dim(x))[2]
                                                int[k,l] <- z[ind.x, ind.y]
                                          } else {
                                                ind.x <- which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[l]) ^ 2))
                                                int[k,l] <- z[ind.x]
                                          }
                                    }
                              }
                        }
                  }
                  z <- NULL
                  return(int)
            })
            if (any(attr(new.Coordinates,"type") == "location")) {
                  obj$Data <- t(unname(do.call("abind", c(interp.list, along = 2))))
                  attr(obj$Data, "dimensions") <- c("time", "station")        
            } else {
                  obj$Data <- unname(do.call("abind", c(interp.list, along = 3)))
                  attr(obj$Data, "dimensions") <- c("lon", "lat", "time")
            }
            interp.list <- NULL
      } else {
            message("[", Sys.time(), "] Performing multi-member ", method, " interpolation... may take a while")
            aux.list <- list()
            # To reduce the copy & paste code we use a kind of wrapper
            # function for lapply 
            if (parallel.pars$hasparallel) {
                apply_fun    <- function(...) {
                    parallel::parLapply(cl = parallel.pars$cl, ...)
                }  
                on.exit(parallel::stopCluster(parallel.pars$cl))
            } else {
                apply_fun <- lapply
            }
            for (i in 1:n.members) {
                  message("[", Sys.time(), "] Interpolating member ", i, " out of ", n.members)
                  interp.list <- apply_fun(1:dim(obj$Data)[time.ind], function(j) {
                        indices <- rep(list(bquote()), length(dim(obj$Data)))
                        indices[[time.ind]] <- j
                        indices[[mem.ind]] <- i
                        ind.call <- as.call(c(list(as.symbol("["), quote(obj$Data)), indices))
                        z <- eval(ind.call)
                        if (isTRUE(transpose)) {
                              z <- t(z)      
                        }
                        indNoNA <- which(!is.nan(z) | !is.na(z))
                        if (any(attr(new.Coordinates,"type") == "location")) {
                              int <- array(data = NA, dim = length(new.Coordinates$x))
                        }else{
                              int <- matrix(data = NA, nrow = length(new.Coordinates$x), ncol = length(new.Coordinates$y))
                        }
                        if (method == "bilinear" & length(indNoNA) != 0) {
                              int <- interp(x = x[indNoNA], y = y[indNoNA], z[indNoNA], xo = new.Coordinates$x, yo = new.Coordinates$y,
                                            linear = TRUE, extrap = FALSE, duplicate = "error",
                                            dupfun = NULL, ncp = NULL, nx = length(new.Coordinates$x), ny = length(new.Coordinates$y))
                              if (any(attr(new.Coordinates,"type") == "location")) {
                                    int <- int[c(0:(length(new.Coordinates$y) - 1)) * length(new.Coordinates$y) + c(1:length(new.Coordinates$y))]
                              }
                              int = int$z
                        }
                        if (method == "nearest" & length(indNoNA) != 0) {
                              if (any(attr(new.Coordinates,"type") == "location")) {
                                    int <- array(data = NA, dim = length(new.Coordinates$x))
                                    for (k in 1:length(new.Coordinates$x)) {
                                          if (!any(attr(obj$Data, "dimensions") == "station")) {
                                                ind.x <- arrayInd(which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[k]) ^ 2)), dim(x))[1]
                                                ind.y <- arrayInd(which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[k]) ^ 2)), dim(x))[2]
                                                int[k] <- z[ind.x, ind.y]
                                          } else {
                                                ind.x <- which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[k]) ^ 2))
                                                int[k] <- z[ind.x]
                                          }
                                    }
                              } else {
                                    int <- matrix(nrow = length(new.Coordinates$x), ncol = length(new.Coordinates$y))
                                    for (k in 1:length(new.Coordinates$x)) {
                                          for (l in 1:length(new.Coordinates$y)) {
                                                if (!any(attr(obj$Data, "dimensions") == "station")) {
                                                      ind.x <- arrayInd(which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[l]) ^ 2)), dim(x))[1]
                                                      ind.y <- arrayInd(which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[l]) ^ 2)), dim(x))[2]
                                                      int[k,l] <- z[ind.x, ind.y]
                                                } else {
                                                      ind.x <- which.min(sqrt((x - new.Coordinates$x[k]) ^ 2 + (y - new.Coordinates$y[l]) ^ 2))
                                                      int[k,l] <- z[ind.x]
                                                }
                                          }
                                    }
                              }
                        }
                        z <- NULL
                        return(int)
                  })

                  if (any(attr(new.Coordinates,"type") == "location")) {
                        aux.list[[i]] <- unname(do.call("abind", c(interp.list, along = 2)))
                  } else {
                        aux.list[[i]] <- unname(do.call("abind", c(interp.list, along = 3)))
                  }
                  interp.list <- NULL
            }
            if (any(attr(new.Coordinates,"type") == "location")) {
                  obj$Data <- unname(do.call("abind", c(aux.list, along = 3)))
                  attr(obj$Data, "dimensions") <- c("station", "time", "member")
            }else{
                  obj$Data <- unname(do.call("abind", c(aux.list, along = 4)))
                  attr(obj$Data, "dimensions") <- c("lon", "lat", "time", "member")
            }
            aux.list <- NULL
      }      
      # Dimension ordering & Coordinate system
      if (!any(attr(new.Coordinates,"type") == "location")) {
            tab <- c("member", "time", "level", "lat", "lon")
            obj$xyCoords$x <- new.Coordinates$x
            obj$xyCoords$y <- new.Coordinates$y
            attr(obj$xyCoords, "resX") <- abs(new.Coordinates$x[2] - new.Coordinates$x[1])
            attr(obj$xyCoords, "resY") <- abs(new.Coordinates$y[2] - new.Coordinates$y[1]) 
            if (is.null(attr(obj$xyCoords, "projection")) & !is.null(attr(new.Coordinates, "projection"))){
                  attr(obj$xyCoords, "projection") <- attr(new.Coordinates, "projection")
            }
      } else {
            tab <- c("member", "time", "level", "station")
            obj$xyCoords <- matrix(abind(new.Coordinates$x,new.Coordinates$y), nrow = length(new.Coordinates$x), ncol = 2)
            attr(obj$xyCoords, "type") <- "location"
            if (!is.null(attr(new.Coordinates, "projection"))){
                  attr(obj$xyCoords, "projection") <- attr(new.Coordinates, "projection")
            }
            
      }
      attr(obj$xyCoords, "interpolation") <-  method
      x <- attr(obj$Data, "dimensions")
      b <- na.exclude(match(tab, x))
      x <- x[b]
      obj$Data <- aperm(obj$Data, perm = b)    
      attr(obj$Data, "dimensions")  <- x
      message("[", Sys.time(), "] Done")
      return(obj)
}
# End

