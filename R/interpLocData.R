#' @title Interpolate a dataset to a grid
#' 
#' @description Performs interpolation of a gridded/location dataset into a new user-defined location using bilinear weights 
#' or nearest-neighbour methods.
#' 
#' @importFrom fields interp.surface
#' @importFrom abind abind
#' 
#' @param obj A data object coming from \code{\link{loadGridData}}, \code{\link{loadStationData}} or the \pkg{ecomsUDG.Raccess} 
#' package function \code{\link[ecomsUDG.Raccess]{loadECOMS}}.
#' @param new.points Definition of the new locations, in the form of a list with the x and y components, in thir order.
#' Each component consists of a vector of length three with components \emph{from}, \emph{to} and \emph{by},
#'  in this order, similar as the arguments passed to the \code{\link[base]{seq}} function, giving the 
#'  westernmost, easternmost and grid cell width in the X axis and, in the same way,
#'  the southernmost, northernmost and grid cell resolution in the Y axis. See details.
#' @param method Method for interpolation. Currently implemented methods are either \code{bilinear},
#' for bilinear interpolation, and \code{nearest}, for nearest-neighbor interpolation (default).
#' @return An interpolated object preserving the output structure of the input
#'  (See e.g. \code{\link{loadGridData}} for details on the output structure). 
#' @details  In case of default definition of either x, y or both coordinates, the default location
#' is the same than the locations considered in the data object \code{obj}.
#' The bilinear interpolator is a wrapper of the \code{\link[fields]{interp.surface}} function
#' in package \pkg{fields}.
#' The output has special attributes in the \code{xyCoords} element that indicate that the object
#'  has been interpolated. These attributes are \code{interpolation}, which indicates the method used and
#'  \code{resX} and \code{resY}, for the grid-cell resolutions in the X and Y axes respectively.
#'  It is also possible to pass the interpolator the grid of a previously existing grid dataset using the
#'  \code{\link{getGrid}} method.
#' @note To avoid unnecessary NA values, the function will not extrapolate using a new grid outside the
#' current extent of the dataset, returning an error message.
#' @family loading.grid
#' @author J. Bedia \email{joaquin.bedia@@gmail.com} and S. Herrera

interpLocData <- function(obj, new.points = list(x = NULL, y = NULL), method = c("nearest", "bilinear")) {
  method <- match.arg(method, choices = c("nearest", "bilinear"))
  if (any(attr(obj$Data, "dimensions") == "station")){
    x <- as.numeric(obj$xyCoords[,1])
    y <- as.numeric(obj$xyCoords[,2])
  }else{
    x <- obj$xyCoords$x
    y <- obj$xyCoords$y
  }
  if (is.null(new.points)) {
    new.points <- getGrid(obj)
  }else{
    if (is.null(new.points$x)) {
      new.points$x <- x
    }
    if (is.null(new.points$y)) {
      new.points$y <- y
    }
  }
  # Definition of new points
  if ((max(new.points$x) < min(x)) | (min(new.points$x) > max(x))) {
    stop("The input and output grids do not overlap\nCheck the input and output grid definitions")
  }
  if (min(new.points$x) < floor(min(x)) | max(new.points$x) > ceiling(max(x))) {
    warning("The new longitudes are outside the data extent")
  }
  if ((max(new.points$y) < min(y)) | (min(new.points$y) > max(y))) {
    stop("The input and output grids do not overlap\nCheck the input and output grid definitions")
  }
  if (min(new.points$y) < floor(min(y)) | max(new.points$y) > ceiling(max(y))) {
    warning("The new latitudes are outside the data extent")
  }
  loc <- matrix(abind(new.points$x,new.points$y), nrow =length(new.points$x), ncol = 2)
  if (any(grepl("member", attr(obj$Data, "dimensions")))) {
    mem.ind <- grep("member", attr(obj$Data, "dimensions"))
    n.members <- dim(obj$Data)[mem.ind]
  } else {
    n.members <- NULL
  }
  time.ind <- grep("^time", attr(obj$Data, "dimensions"))
  if (!any(attr(obj$Data, "dimensions") == "station")){
    # Handles reverse ordering of lon-lat (i.e. lat-lon)
    ind.order <- match(c("lon", "lat"), attr(obj$Data, "dimensions"))
    transpose <- ifelse(!identical(sort(ind.order), ind.order), TRUE, FALSE)
  }else{
    transpose <- FALSE
  }
  if (is.null(n.members)) {
    message("[", Sys.time(), "] Performing ", method, " interpolation... may take a while")
    interp.list <- lapply(1:dim(obj$Data)[time.ind], function(j) {
      indices <- rep(list(bquote()), length(dim(obj$Data)))
      indices[[time.ind]] <- j
      ind.call <- as.call(c(list(as.symbol("["), quote(obj$Data)), indices))
      z <- eval(ind.call)
      if (isTRUE(transpose)) {
        z <- t(z)      
      }
      if (method == "bilinear") {
        int <- interp.surface(list("x" = x, "y" = y, "z" = z), loc)
      }
      if (method == "nearest") {
        int <- array(data = NA, dim = length(new.points$x))
        for (k in 1:length(new.points$x)) {
          if (!any(attr(obj$Data, "dimensions") == "station")){
            ind.x <- which.min(abs(x - loc[k,1]))
            ind.y <- which.min(abs(y - loc[k,2]))
            int[k] <- z[ind.x, ind.y]
          }else{
            ind.x <- which.min(sqrt((x - loc[k,1])^2+(y - loc[k,2])^2))
            int[k] <- z[ind.x]
          }
        }
      }
      z <- NULL
      return(int)
    })
    obj$Data <- t(unname(do.call("abind", c(interp.list, along = 2))))
    interp.list <- NULL
    attr(obj$Data, "dimensions") <- c("time", "station")
  } else {
    message("[", Sys.time(), "] Performing multi-member ", method, " interpolation... may take a while")
    aux.list <- list()
    for (i in 1:n.members) {
      message("[", Sys.time(), "] Interpolating member ", i, " out of ", n.members)
      interp.list <- lapply(1:dim(obj$Data)[time.ind], function(j) {
        indices <- rep(list(bquote()), length(dim(obj$Data)))
        indices[[time.ind]] <- j
        indices[[mem.ind]] <- i
        ind.call <- as.call(c(list(as.symbol("["), quote(obj$Data)), indices))
        z <- eval(ind.call)
        if (isTRUE(transpose)) {
          z <- t(z)      
        }
        if (method == "bilinear") {
          int <- interp.surface(list("x" = x, "y" = y, "z" = z), loc)
        }
        if (method == "nearest") {
          int <- array(data = NA, dim = length(new.points$x))
          for (k in 1:length(new.points$x)) {
            if (!any(attr(obj$Data, "dimensions") == "station")){
              ind.x <- which.min(abs(x - loc[k,1]))
              ind.y <- which.min(abs(y - loc[k,2]))
              int[k] <- z[ind.x, ind.y]
            }else{
              ind.x <- which.min(sqrt((x - loc[k,1])^2+(y - loc[k,2])^2))
              int[k] <- z[ind.x]
            }
          }
        }
        z <- NULL
        return(int)
      })
      aux.list[[i]] <- unname(do.call("abind", c(interp.list, along = 2)))
    }
    interp.list <- NULL
    obj$Data <- unname(do.call("abind", c(aux.list, along = 3)))
    aux.list <- NULL
    attr(obj$Data, "dimensions") <- c("station", "time", "member")
  }      
  obj$xyCoords <- loc
  attr(obj$xyCoords, "interpolation") <-  method
  attr(obj$xyCoords, "type") <- "location"
  # Dimension ordering
  tab <- c("member", "time", "level", "station")
  x <- attr(obj$Data, "dimensions")
  b <- na.exclude(match(tab, x))
  x <- x[b]
  obj$Data <- aperm(obj$Data, perm = b)    
  attr(obj$Data, "dimensions")  <- x
  message("[", Sys.time(), "] Done")
  return(obj)
}
# End   
