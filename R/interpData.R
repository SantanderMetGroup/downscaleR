#' @title Interpolate a dataset
#' @description Interpolation of gridded datasets into a user-defined grid using nearest-neighbour or bilinear weights. 
#' @importFrom akima interp
#' @importFrom abind abind
#' @importFrom fields interp.surface.grid
#' @param obj A data object coming from \code{loadGridData}, \code{loadStationData} of package \pkg{loadeR}, or the \pkg{loadeR.ECOMS} 
#' package function \code{loadECOMS}.
#' @param new.coordinates Definition of the new coordinates (grid or locations), in the form of a list with the x and y components, in thir order.
#' Each component consists of a vector of length three with components \emph{from}, \emph{to} and \emph{by},
#'  in this order, similar as the arguments passed to the \code{\link[base]{seq}} function, giving the 
#'  westernmost, easternmost and grid cell width in the X axis and, in the same way,
#'  the southernmost, northernmost and grid cell resolution in the Y axis. See details.
#' @param method Method for interpolation. Currently implemented methods are either \code{"bilinear"},
#' for bilinear interpolation, and \code{"nearest"}, for nearest-neighbor interpolation (default).
#' @template templateParallelParams 
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
#' @template templateParallel
#' @note To avoid unnecessary NA values, the function will not extrapolate using a new grid outside the
#' current extent of the dataset, returning an error message.
#' @author J. Bedia, S. Herrera, M. de Felice
#' @export
#' @examples \dontrun{
#' # Load air temperature at 850 mb isobaric pressure level for boreal winter (DJF) 1991-2010
#' data(iberia_ncep_ta850)
#' par(mfrow = c(2,1))
#' plotMeanGrid(iberia_ncep_ta850)
#' # Bilinear interpolation to domain centered in Spain using a 0.5 degree resolution 
#' # in both X and Y axes
#' t <- interpData(iberia_ncep_ta850, new.coordinates = list(x = c(-10,5,.5),
#'                                                           y = c(36,44,.5)),
#'                                    method = "bilinear")
#' plotMeanGrid(t)
#' par(mfrow=c(1,1))
#' # New attributes indicate that the data have been interpolated:
#' attributes(t$xyCoords)
#' }

interpData <- function(obj,
                       new.coordinates = list(x = NULL, y = NULL),
                       method = c("nearest", "bilinear"),
                       parallel = FALSE,
                       max.ncores = 16,
                       ncores = NULL) {
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
                  x <- t(list(x = outer(obj$xyCoords$y*0,obj$xyCoords$x, FUN = "+"),
                              y = outer(obj$xyCoords$y,obj$xyCoords$x*0, FUN = "+"))$x)
                  y <- t(list(x = outer(obj$xyCoords$y*0,obj$xyCoords$x, FUN = "+"),
                              y = outer(obj$xyCoords$y,obj$xyCoords$x*0, FUN = "+"))$y)
            }
      }
      if (is.null(new.coordinates)) {
            new.coordinates <- getGrid(obj)
      } else {
            if (is.null(new.coordinates$x)) {
                  new.coordinates$x <- x
            } else if (exists("resX", where = attributes(new.coordinates))) {
                  if (length(new.coordinates$x) != 2 | new.coordinates$x[2] < new.coordinates$x[1]) {
                        stop("Invalid grid definition in X")
                  }
                  if ((max(c(new.coordinates$x[1],new.coordinates$x[2])) < min(x)) | (min(c(new.coordinates$x[1],new.coordinates$x[2])) > max(x))) {
                        stop("The input and output grids do not overlap\nCheck the input and output grid definitions")
                  }
                  if (new.coordinates$x[1] < floor(min(x)) | new.coordinates$x[2] > ceiling(max(x))) {
                        warning("The new longitudes are outside the data extent")
                  }
                  new.coordinates$x <- do.call("seq", as.list(c(new.coordinates$x, attr(new.coordinates, 'resX'))))
            } else if (length(new.coordinates$x) == 3) {
                  if (new.coordinates$x[2] < new.coordinates$x[1]) {
                        stop("Invalid grid definition in X")
                  }
                  if ((max(c(new.coordinates$x[1],new.coordinates$x[2])) < min(x)) | (min(c(new.coordinates$x[1],new.coordinates$x[2])) > max(x))) {
                        stop("The input and output grids do not overlap\nCheck the input and output grid definitions")
                  }
                  if (new.coordinates$x[1] < floor(min(x)) | new.coordinates$x[2] > ceiling(max(x))) {
                        warning("The new longitudes are outside the data extent")
                  }
                  if ((new.coordinates$x[2] > new.coordinates$x[1]) & (abs(new.coordinates$x[3]) < abs(new.coordinates$x[2] - new.coordinates$x[1]))) {
                        new.coordinates$x <- seq(from = new.coordinates$x[1], to = new.coordinates$x[2], by = new.coordinates$x[3])
                  }
            }
            if (is.null(new.coordinates$y)) {
                  new.coordinates$y <- y
            } else if (exists("resY", where = attributes(new.coordinates))) {
                  if (length(new.coordinates$y) != 2 | new.coordinates$y[2] < new.coordinates$y[1]) {
                        stop("Invalid grid definition in Y")
                  }
                  if ((max(c(new.coordinates$y[1],new.coordinates$y[2])) < min(y)) | (min(c(new.coordinates$y[1],new.coordinates$y[2])) > max(y))) {
                        stop("The input and output grids do not overlap\nCheck the input and output grid definitions")
                  }
                  if (new.coordinates$y[1] < floor(min(y)) | new.coordinates$y[2] > ceiling(max(y))) {
                        warning("The new latitudes are outside the data extent")
                  }
                  new.coordinates$y <- do.call("seq", as.list(c(new.coordinates$y, attr(new.coordinates, 'resY'))))
            } else if (length(new.coordinates$y) == 3) {
                  if (new.coordinates$y[2] < new.coordinates$y[1]) {
                        stop("Invalid grid definition in Y")
                  }
                  if ((max(c(new.coordinates$y[1],new.coordinates$y[2])) < min(y)) | (min(c(new.coordinates$y[1],new.coordinates$y[2])) > max(y))) {
                        stop("The input and output grids do not overlap\nCheck the input and output grid definitions")
                  }
                  if (new.coordinates$y[1] < floor(min(y)) | new.coordinates$y[2] > ceiling(max(y))) {
                        warning("The new latitudes are outside the data extent")
                  }
                  if ((new.coordinates$y[2] > new.coordinates$y[1]) & (abs(new.coordinates$y[3]) < abs(new.coordinates$y[2] - new.coordinates$y[1]))) {
                        new.coordinates$y <- seq(from = new.coordinates$y[1], to = new.coordinates$y[2], by = new.coordinates$y[3])
                  }
            }
      }
      # function for lapply 
      if (parallel.pars$hasparallel) {
            apply_fun <- function(...) {
                  parallel::parLapply(cl = parallel.pars$cl, ...)
            }  
            on.exit(parallel::stopCluster(parallel.pars$cl))
      } else {
            apply_fun <- lapply
      }
      # redim object
      obj <- redim(obj, runtime = FALSE)
      mem.ind <- grep("member", getDim(obj))
      n.members <- dim(obj$Data)[mem.ind]
      time.ind <- grep("^time", getDim(obj))
      n.times <- dim(obj$Data)[time.ind]
      # nearest indices
      if (method == "nearest") {
            if (any(attr(new.coordinates, "type") == "location")) {
                  ind.NN <- array(data = NA, dim = c(1,length(new.coordinates$x)))
                  for (k in 1:length(new.coordinates$x)) {
                        distK <- sqrt((x - new.coordinates$x[k]) ^ 2 + (y - new.coordinates$y[k]) ^ 2)
                        ind.NN[k] <- which.min(distK)
                  }
            } else {
                  ind.NN <- matrix(nrow = length(new.coordinates$x), ncol = length(new.coordinates$y))
                  for (k in 1:length(new.coordinates$x)) {
                        for (l in 1:length(new.coordinates$y)) {
                              distK <- sqrt((x - new.coordinates$x[k]) ^ 2 + (y - new.coordinates$y[l]) ^ 2)
                              ind.NN[k,l] <- which.min(distK)
                        }
                  }
            }
      }
      message("[", Sys.time(), "] Performing ", method, " interpolation... may take a while")
      aux.list <- list()
      for (i in 1:n.members) {
            if (n.members > 1) message("[", Sys.time(), "] Interpolating member ", i, " out of ", n.members)
            if (method == "nearest") {
                  if (any(attr(new.coordinates,"type") == "location")) {
                        int <- array(data = NA, dim = c(n.members,length(new.coordinates$x)))
                        for (k in 1:length(new.coordinates$x)) {
                              if (!any(attr(obj$Data, "dimensions") == "station")) {
                                    ind.x <- arrayInd(ind.NN[k], dim(x))[1]
                                    ind.y <- arrayInd(ind.NN[k], dim(x))[2]
                                    int[,k] <- obj$Data[,ind.y, ind.x]
                              } else {
                                    int[,k] <- obj$Data[,ind.NN[k]]
                              }
                        }
                  } else {
                        int <- array(dim = c(n.times, length(new.coordinates$y), length(new.coordinates$x)))
                        for (k in 1:length(new.coordinates$x)) {
                              for (l in 1:length(new.coordinates$y)) {
                                    if (!any(attr(obj$Data, "dimensions") == "station")) {
                                          ind.x <- arrayInd(ind.NN[k,l], dim(x))[1]
                                          ind.y <- arrayInd(ind.NN[k,l], dim(x))[2]
                                          int[,l,k] <- obj$Data[i,,ind.y, ind.x]
                                    } else {
                                          int[,l,k] <- obj$Data[,,ind.NN[k,l]]
                                    }
                              }
                        }
                  }
                  aux.list[[i]] <- int
                  int <- NULL
                  dimNames.ref <- c("member","time","lat","lon")
            }
            if (method == "bilinear") {
                  interp.list <- apply_fun(1:n.times, function(j) { # iterates in time (inefficient!, to be changed)
                        z <- asub(obj$Data, idx = list(i,j), dims = c(mem.ind, time.ind))
                        any_is_NA_or_NAN <- any(is.nan(z) | anyNA(z))
                        if (any_is_NA_or_NAN) message("NOTE: The presence of NAs/NANs may slow down interpolation")
                        int <- if (any(attr(new.coordinates,"type") == "location")) {
                              array(data = NA, dim = length(new.coordinates$x))
                        } else {
                              matrix(nrow = length(new.coordinates$x), ncol = length(new.coordinates$y))
                        }
                        if (any_is_NA_or_NAN) {
                              # Due to the presence of NAs or NANs akima::interp must be used (slower)
                              indNoNA <- which(!is.nan(z) | !is.na(z))
                              int <- interp(x = x[indNoNA], y = y[indNoNA], z[indNoNA], xo = new.coordinates$x, yo = new.coordinates$y,
                                            linear = TRUE, extrap = FALSE, duplicate = "error",
                                            dupfun = NULL, ncp = NULL, nx = length(new.coordinates$x), ny = length(new.coordinates$y))
                              if (any(attr(new.coordinates,"type") == "location")) {
                                    int <- int[c(0:(length(new.coordinates$y) - 1)) * length(new.coordinates$y) + c(1:length(new.coordinates$y))]
                              } else {
                                    int = int$z    
                              }
                        } else { 
                              # In absence of NAs/NANs we use fields::interp.surface.grid faster than akima::interp
                              int <- fields::interp.surface.grid(list(x = obj$xyCoords$x,
                                                                      y = obj$xyCoords$y,
                                                                      z = t(z)),
                                                                 grid.list = list(x = new.coordinates$x,
                                                                                  y = new.coordinates$y))
                              int <- if (any(attr(new.coordinates,"type") == "location")) {
                                    int[c(0:(length(new.coordinates$y) - 1)) * length(new.coordinates$y) + c(1:length(new.coordinates$y))]
                              } else {
                                    int$z    
                              }
                        }
                        z <- NULL
                        return(int)
                  })
                  if (any(attr(new.coordinates,"type") == "location")) {
                        aux.list[[i]] <- unname(do.call("abind", c(interp.list, along = 2)))
                  } else {
                        aux.list[[i]] <- unname(do.call("abind", c(interp.list, along = -1L)))
                  }
                  interp.list <- NULL
                  dimNames.ref <- c("member", "time", "lon", "lat")
            }
      }
      if (any(attr(new.coordinates,"type") == "location")) {
            obj$Data <- unname(do.call("abind", c(aux.list, along = 3L)))
            attr(obj$Data, "dimensions") <- c("station", "time", "member")
      } else {
            obj$Data <- unname(do.call("abind", c(aux.list, along = -1L)))
            attr(obj$Data, "dimensions") <- dimNames.ref
            dimNames.ref <- c("member", "time", "lon", "lat")
            }
      aux.list <- NULL
      # Dimension ordering & Coordinate system
      if (!any(attr(new.coordinates,"type") == "location")) {
            tab <- c("member", "time", "level", "lat", "lon")
            obj$xyCoords$x <- new.coordinates$x
            obj$xyCoords$y <- new.coordinates$y
            attr(obj$xyCoords, "resX") <- abs(new.coordinates$x[2] - new.coordinates$x[1])
            attr(obj$xyCoords, "resY") <- abs(new.coordinates$y[2] - new.coordinates$y[1]) 
            if (is.null(attr(obj$xyCoords, "projection")) & !is.null(attr(new.coordinates, "projection"))) {
                  attr(obj$xyCoords, "projection") <- attr(new.coordinates, "projection")
            }
      } else {
            tab <- c("member", "time", "level", "station")
            obj$xyCoords <- matrix(abind(new.coordinates$x,new.coordinates$y), nrow = length(new.coordinates$x), ncol = 2)
            attr(obj$xyCoords, "type") <- "location"
      }
      attr(obj$xyCoords, "interpolation") <-  method
      # x <- getDim(obj)
      # b <- na.exclude(match(tab, x))
      # x <- x[b]
      # obj$Data <- aperm(obj$Data, perm = )
      # attr(obj$Data, "dimensions")  <- x
      if (is.null(attr(obj$xyCoords, "projection"))) {
            attr(obj$xyCoords, "projection") <- "undefined"
      }
      obj <- redim(obj, drop = TRUE)
      message("[", Sys.time(), "] Done")
      return(obj)
}
# End


