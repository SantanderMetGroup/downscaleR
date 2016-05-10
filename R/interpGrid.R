#' @title Grid interpolation
#' @description Interpolation of gridded datasets into a user-defined grid using nearest-neighbour or bilinear weights. 
#' @importFrom akima interp
#' @importFrom abind abind
#' @importFrom fields interp.surface.grid
#' @param grid An input grid to be interpolated/regridded.
#' @param new.coordinates Definition of the new grid coordinates, in the form of a list with the x and y components, in thir order.
#' Each component consists of a vector of length three with components \emph{from}, \emph{to} and \emph{by},
#'  in this order, similar as the arguments passed to the \code{\link[base]{seq}} function, giving the 
#'  westernmost, easternmost and grid cell width in the X axis and, in the same way,
#'  the southernmost, northernmost and grid cell resolution in the Y axis. See details.
#' @param method Method for interpolation. Currently implemented methods are either \code{"bilinear"},
#' for bilinear interpolation, and \code{"nearest"}, for nearest-neighbor interpolation (default).
#' @param bilin.method Algorithm chosen for bilinear interpolation. Two options available: \code{"akima"} uses \code{\link[akima]{interp}} and
#' \code{"fields"} the \code{\link[fields]{interp.surface.grid}} algorithm. In case any missing values exist in the input data matrix, 
#' the \code{"fields"} option, able to handle missing values, need to be used. Otherwise, the \code{"akima"} option performs much faster.
#' @template templateParallelParams 
#' @return An interpolated object preserving the structure of the input
#' @details  In case of default definition of either x, y or both grid coordinates, the default grid
#' is calculated taking the corners of the current grid and assuming x and y resolutions equal to 
#' the default \code{by} argument value in function \code{\link[base]{seq}}: \emph{by = ((to - from)/(length.out - 1))}.
#' The output has special attributes in the \code{xyCoords} element that indicate that the object
#'  has been interpolated. These attributes are \code{interpolation}, which indicates the method used and
#'  \code{resX} and \code{resY}, for the grid-cell resolutions in the X and Y axes respectively.
#'  It is also possible to pass the interpolator the grid of a previously existing grid dataset using the
#'  \code{\link{getGrid}} method.
#' @template templateParallel
#' @note To avoid unnecessary NA values, the function will not extrapolate using a new grid outside the
#' current extent of the dataset, returning an error message.
#' @author J. Bedia, S. Herrera, M. de Felice, M. Iturbide
#' @export
#' @examples \dontrun{
#' # Load air temperature at 850 mb isobaric pressure level for boreal winter (DJF) 1991-2010
#' data(iberia_ncep_ta850)
#' par(mfrow = c(2,1))
#' plotMeanGrid(iberia_ncep_ta850)
#' # Bilinear interpolation to domain centered in Spain using a 0.5 degree resolution 
#' # in both X and Y axes
#' t <- interpGrid(iberia_ncep_ta850, new.coordinates = list(x = c(-10,5,.5),
#'                                                           y = c(36,44,.5)),
#'                                    method = "bilinear",
#'                                    bilin.method = "akima")
#' plotMeanGrid(t)
#' par(mfrow=c(1,1))
#' # New attributes indicate that the data have been interpolated:
#' attributes(t$xyCoords)
#' }

interpGrid <- function(grid,
                       new.coordinates = list(x = NULL, y = NULL),
                       method = c("nearest", "bilinear"),
                       bilin.method = NULL,
                       parallel = FALSE,
                       max.ncores = 16,
                       ncores = NULL) {
      method <- match.arg(method, choices = c("nearest", "bilinear"))
      if (method == "nearest" & !is.null(bilin.method)) message("NOTE: argument 'bilin.method' ignored for nearest neighbour interpolation")
      if (method == "bilinear") bilin.method <- match.arg(bilin.method, choices = c("akima", "fields"))
      parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
      x <- grid$xyCoords$x
      y <- grid$xyCoords$y
      if (!is.null(attr(grid$xyCoords, which = "projection"))) {
            if (attr(grid$xyCoords, which = "projection") == "RotatedPole") {
                  x <- grid$xyCoords$lon
                  y <- grid$xyCoords$lat
            }
      }
      if (is.vector(x) & is.vector(y)) {
            x <- t(list(x = outer(grid$xyCoords$y*0,grid$xyCoords$x, FUN = "+"),
                        y = outer(grid$xyCoords$y,grid$xyCoords$x*0, FUN = "+"))$x)
            y <- t(list(x = outer(grid$xyCoords$y*0,grid$xyCoords$x, FUN = "+"),
                        y = outer(grid$xyCoords$y,grid$xyCoords$x*0, FUN = "+"))$y)
      }
      if (is.null(new.coordinates)) {
            new.coordinates <- getGrid(grid)
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
      grid <- redim(grid, runtime = FALSE)
      mem.ind <- grep("member", getDim(grid))
      n.members <- dim(grid$Data)[mem.ind]
      time.ind <- grep("^time", getDim(grid))
      n.times <- dim(grid$Data)[time.ind]
      # nearest indices
      if (method == "nearest") {
            ind.NN <- matrix(nrow = length(new.coordinates$x), ncol = length(new.coordinates$y))
            for (k in 1:length(new.coordinates$x)) {
                  for (l in 1:length(new.coordinates$y)) {
                        distK <- sqrt((x - new.coordinates$x[k]) ^ 2 + (y - new.coordinates$y[l]) ^ 2)
                        ind.NN[k,l] <- which.min(distK)
                  }
            }
      }
      message("[", Sys.time(), "] Performing ", method, " interpolation... may take a while")
      aux.list <- list()
      for (i in 1:n.members) {
            if (n.members > 1) message("[", Sys.time(), "] Interpolating member ", i, " out of ", n.members)
            if (method == "nearest") {
                  int <- array(dim = c(n.times, length(new.coordinates$y), length(new.coordinates$x)))
                  for (k in 1:length(new.coordinates$x)) {
                        for (l in 1:length(new.coordinates$y)) {
                              ind.x <- arrayInd(ind.NN[k,l], dim(x))[1]
                              ind.y <- arrayInd(ind.NN[k,l], dim(x))[2]
                              int[,l,k] <- grid$Data[i,,ind.y, ind.x]
                        }
                  }
                  aux.list[[i]] <- int
                  int <- NULL
                  dimNames.ref <- c("member","time","lat","lon")
            }
            if (method == "bilinear") {
                  dimNames.ref <- c("member", "time", "lon", "lat")
                  interp.list <- apply_fun(1:n.times, function(j) { # iterates in time (inefficient!, to be changed)
                        z <- asub(grid$Data, idx = list(i,j), dims = c(mem.ind, time.ind))
                        any_is_NA_or_NAN <- any(!is.finite(z))
                        if (bilin.method == "akima") {
                              if (any_is_NA_or_NAN) stop("The input grid contains missing values\nConsider using 'bilin.method=\"fields\"' instead", call. = FALSE)
                              indNoNA <- which(is.finite(z))
                              int <- akima::interp(x = x[indNoNA], y = y[indNoNA], t(z)[indNoNA],
                                                   xo = new.coordinates$x, yo = new.coordinates$y,
                                                   linear = TRUE, extrap = FALSE, duplicate = "error",
                                                   nx = length(new.coordinates$x), ny = length(new.coordinates$y))$z
                        } else if (bilin.method == "fields") {
                              if (!any_is_NA_or_NAN & i == 1 & j == 1) message("NOTE: No missing values present in the input grid\nConsider using 'bilin.method=\"akima\"' for improved speed")
                              int <- fields::interp.surface.grid(list(x = grid$xyCoords$x,
                                                                      y = grid$xyCoords$y,
                                                                      z = t(z)),
                                                                 grid.list = list(x = new.coordinates$x,
                                                                                  y = new.coordinates$y))$z
                        }
                        z <- NULL
                        return(int)
                  })
                  aux.list[[i]] <- unname(do.call("abind", c(interp.list, along = -1L)))
                  interp.list <- NULL
            }
      }
      grid$Data <- unname(do.call("abind", c(aux.list, along = -1L)))
      attr(grid$Data, "dimensions") <- dimNames.ref
      aux.list <- NULL
      # Dimension ordering & Coordinate system
      tab <- c("member", "time", "level", "lat", "lon")
      grid$xyCoords$x <- new.coordinates$x
      grid$xyCoords$y <- new.coordinates$y
      attr(grid$xyCoords, "resX") <- abs(new.coordinates$x[2] - new.coordinates$x[1])
      attr(grid$xyCoords, "resY") <- abs(new.coordinates$y[2] - new.coordinates$y[1]) 
      if (is.null(attr(grid$xyCoords, "projection")) & !is.null(attr(new.coordinates, "projection"))) {
            attr(grid$xyCoords, "projection") <- attr(new.coordinates, "projection")
      }
      attr(grid$xyCoords, "interpolation") <-  method
      x <- getDim(grid)
      b <- na.exclude(match(tab, x))
      x <- x[b]
      grid$Data <- aperm(grid$Data, perm = b)
      attr(grid$Data, "dimensions")  <- x
      if (is.null(attr(grid$xyCoords, "projection"))) {
            attr(grid$xyCoords, "projection") <- "undefined"
      }
      grid <- redim(grid, drop = TRUE)
      message("[", Sys.time(), "] Done")
      return(grid)
}
# End


