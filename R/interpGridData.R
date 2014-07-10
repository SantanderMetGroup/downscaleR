#' Interpolate a dataset to a grid
#' 
#' Usaes bilinear weights or nearest neighbour interpolation for interpolating a gridded
#'  dataset into a new grid.
#'  
#'  @param obj An object coming from \code{\link{loadGridData}} or the \code{ecomsUDG.Raccess} package function
#'   \code{\link[ecomsUDG.Raccess]{loadECOMS}}.
#'  @param new.grid.x Definition of the x coordinates of the grid to interpolate.
#'  This is a vector of length three with components \emph{from}, \emph{to} and \emph{by},
#'   in this order, similar as the arguments passed to the \code{\link[base]{seq}} function, giving the 
#'   westernmost, easternmost and grid cell width in the X axis parameters. See details.
#'  @param new.grid.y Same as \code{new.grid.x} but for the Y coordinates, giving the southernmost,
#'   northernmost and grid cell resolution in the Y axis. See details
#'  @param method Currently two interpolation methods are accepted: bilinear interpolation (default)
#'   and nearest neighbour (faster for large data).
#'  @return An interpolated object preserving the output structure of the input
#'  @details  In case of default definition of either x, y or both grid coordinates, the default grid
#'  is calculated taking the corners of the current grid and assuming x and y resolutions equal to 
#'  the default \code{by} arguent value in function \code{\link[base]{seq}}: \emph{by = ((to - from)/(length.out - 1))}.
#'  The bilinear interpolator is essentially a wrapper of the \code{fields} package 
#'  function \code{\link[fields]{interp.surface.grid}}.
#'  @note To avoid unnecessary NA values within the dataset, the function will not accept new grid domains outside the
#'  current extent of the dataset, returning an error message.
#'  @author J. Bedia \email{joaquin.bedia@@gmail.com}
#'  @export

# gridData <- loadGridData(dataset = "inst/datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.ncml", 
#                          var = "ta@100000",
#                          dictionary = TRUE,
#                          lonLim = c(-12, 5),
#                          latLim = c(35,45),
#                          season = 6:8,
#                          years = 1981:2000,
#                          time = "none")
# 
# 
# library(ecomsUDG.Raccess)
# loginECOMS_UDG("juaco", "*********")
# 
# gridData <- loadECOMS(dataset = "System4_seasonal_15",
#           var = "psl",
#           dictionary = TRUE,
#           members = 1:4,
#           lonLim = c(-15,15),
#           latLim = c(33,45),
#           season = 6:9,
#           years = 1999:2003,
#           leadMonth = 2,
#           time = "12") 
#             
# 
# plotMeanField(gridData)
# 
# new.grid.x <- c(-10,5,.25)
# new.grid.y <- c(35,43,.25)
# 
# new.grid.x <- c(-16,38,.25)
# 
# x
# y

interpGridData <- function(gridData, new.grid.x = NULL, new.grid.y = NULL, method = c("bilinear", "nearest")) {
      x <- gridData$xyCoords$x
      y <- gridData$xyCoords$y
      # Definition of new grid
      if (is.null(new.grid.x)) {
            new.grid.x <- seq(x[1], tail(x, 1))
      } else {
            if (length(new.grid.x) != 3 | new.grid.x[2] < new.grid.x[1]) {
                  stop("Invalid grid definition in X")
            }
            if (new.grid.x[1] < floor(x[1]) & new.grid.x[2] <= ceiling(tail(x, 1))) {
                  stop("The westernmost corner of the new grid is outside the data extent\n Minimum X accepted value: ", floor(x[1]))
            }
            if (new.grid.x[2] > ceiling(tail(x, 1)) & new.grid.x[1] >= floor(x[1])) {
                  stop("The easternmost corner of the new grid is outside the data extent\n Maximum X accepted value: ", ceiling(tail(x, 1)))
            }
            if (new.grid.x[2] > ceiling(tail(x, 1)) & new.grid.x[1] < floor(x[1])) {
                  stop("The new grid is outside the data extent\n Accepted X values in the range: [", floor(x[1]), ",", ceiling(tail(x, 1)), "]")
            }
            new.grid.x <- do.call("seq", as.list(new.grid.x))
      }
      if (is.null(new.grid.y)) {
            new.grid.y <- seq(y[1], tail(y, 1))
      } else {
            if (length(new.grid.y) != 3 | new.grid.y[2] < new.grid.y[1]) {
                  stop("Invalid grid definition in Y")
            }
            
            if (new.grid.y[1] < floor(y[1]) & new.grid.y[2] <= ceiling(tail(y, 1))) {
                  stop("The southernmost corner of the new grid is outside the data extent\n Minimum Y accepted value: ", floor(y[1]))
            }
            if (new.grid.y[2] > ceiling(tail(y, 1)) & new.grid.y[1] >= floor(y[1])) {
                  stop("The northernmost corner of the new grid is outside the data extent\n Maximum Y accepted value: ", ceiling(tail(y, 1)))
            }
            if (new.grid.y[2] > ceiling(tail(y, 1)) & new.grid.y[1] < floor(y[1])) {
                  stop("The new grid is outside the data extent\n Accepted Y values in the range: [", floor(y[1]), ",", ceiling(tail(y, 1)), "]")
            }
            new.grid.y <- do.call("seq", as.list(new.grid.y))
      }
      grid.list <- list("x" = new.grid.x, "y" = new.grid.y)
      new.grid.x <- NULL
      new.grid.y <- NULL


#       
#       # nearest

# which.min(abs(x-your.number))
#       ng <- expand.grid(grid.list$x, grid.list$y)
#       og <- expand.grid(x, y)
#       
#       for (i in 1:nrow(ng)) {
#                         
#       }
#             
#       
#       for (i in 1:length(grid.list$x))
#             
#             i=100
#             
#             grid.list$x[i]^2
      
      # Bilinear
      if (any(grepl("member", attr(gridData$Data, "dimensions")))) {
            mem.ind <- grep("member", attr(gridData$Data, "dimensions"))
            n.members <- dim(gridData$Data)[mem.ind]
      } else {
            n.members <- NULL
      }
      time.ind <- grep("^time", attr(gridData$Data, "dimensions"))
      # Handles reverse ordering of lon-lat (i.e. lat-lon)
      ind.order <- match(c("lon", "lat"), attr(gridData$Data, "dimensions"))
      transpose <- FALSE
      if (!identical(sort(ind.order), ind.order)) {
            transpose <- TRUE
      }
      if (!is.null(n.members)) {
            message("[", Sys.time(), "] Performing multi-member interpolation (", n.members, " members) ... may take a while")
            aux.list <- list()
            for (i in 1:n.members) {
                  message("[", Sys.time(), "] Interpolating member ", i, " ...")
                  interp.list <- lapply(1:dim(gridData$Data)[time.ind], function(j) {
                        indices <- rep(list(bquote()), length(dim(gridData$Data)))
                        indices[[time.ind]] <- j
                        indices[[mem.ind]] <- i
                        ind.call <- as.call(c(list(as.symbol("["), quote(gridData$Data)), indices))
                        z <- eval(ind.call)
                        if (isTRUE(transpose)) {
                              z <- t(z)      
                        } 
                        obj <- list("x" = x, "y" = y, "z" = z)
                        z <- NULL
                        int <- interp.surface.grid(obj, grid.list)$z
                        obj <- NULL
                        return(int)
                  })
                  aux.list[[i]] <- unname(do.call("abind", c(interp.list, along = 3)))
                  interp.list <- NULL
            }
            gridData$Data <- unname(do.call("abind", c(aux.list, along = 4)))
            aux.list <- NULL
            attr(gridData$Data, "dimensions") <- c("lon", "lat", "time", "member")
      } else {
            message("[", Sys.time(), "] Performing bilinear interpolation... may take a while")
            interp.list <- lapply(1:dim(gridData$Data)[time.ind], function(j) {
                  indices <- rep(list(bquote()), length(dim(gridData$Data)))
                  indices[[time.ind]] <- j
                  indices[[mem.ind]] <- i
                  ind.call <- as.call(c(list(as.symbol("["), quote(gridData$Data)), indices))
                  z <- eval(ind.call)
                  if (isTRUE(transpose)) {
                        z <- t(z)      
                  } 
                  obj <- list("x" = x, "y" = y, "z" = z)
                  z <- NULL
                  int <- interp.surface.grid(obj, grid.list)$z
                  obj <- NULL
                  return(int)
            })
            grid$Data <- unname(do.call("abind", c(interp.list, along = 3)))
            interp.list <- NULL
            attr(gridData$Data, "dimensions") <- c("lon", "lat", "time")
      }      
      message("[", Sys.time(), "] Done")
      gridData$xyCoords[1:2] <- grid.list
      return(gridData)
}
# End   
 

# a <- interpGridData(gridData, new.grid.x = c(-10, 5, .2), new.grid.y = c(35,43,.2), method = "bilinear")
# plotMeanField(a)
# str(a)
            