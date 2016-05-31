#' @title Lattice plot methods for climatological grids
#' @description A wrapper for the lattice (trellis) plot methods for spatial data in \code{sp::spplot}
#' @param grid Input grid
#' @param tolerance Precision, used to which extent points are exactly on a grid. See details.
#' @param backdrop.theme Reference geographical lines to be added to the plot. See Details. 
#' @param ... Further arguments passed to \code{spplot}
#' @details The function applies the \code{\link[sp]{spplot}} method after conversion of the climatological map(s) to a
#'  \code{SpatialGridDataFrame}.
#'  
#'  \strong{Tolerance parameter}
#'   
#'   In case of not perfectly regular grids, an error will raise suggesting a minimum value to be passed to
#'  \code{tolerance}. This argument is passed to \code{\link[sp]{points2grid}}.
#'  
#'  \strong{Backdrop theme}
#'  
#'  Current implemented options are \code{"none"} and \code{"coastline"}, which contains
#'  a simplied vector theme delineating the world coastlines.
#'  
#' @importFrom abind abind
#' @importClassesFrom sp SpatialGridDataFrame SpatialPoints
#' @importFrom sp points2grid gridded coordinates spplot
#' @importFrom grDevices colorRampPalette
#' @export
#' @author J. Bedia

## Lattice plotMeanGrid
# load("~/Desktop/documentation/previous_tools/data/demoData.rda", verbose = TRUE)
# grid <- tx.forecast
# grid <- climatology(tx.forecast, by.member = TRUE)
# grid <- climatology(tx.obs)
# tolerance = 0.004
# 
# 
# save()


plotClimatology <- function(grid,
                            tolerance = sqrt(.Machine$double.eps),
                            backdrop.theme = "none",
                            ...) {
      arg.list <- list(...)
      bt <- match.arg(backdrop.theme, choices = c("none", "coastline", "countries"))
      grid <- redim(grid, drop = FALSE)
      dimNames <- getDim(grid)
      mem.ind <- grep("member", dimNames)
      n.mem <- dim(grid[["Data"]])[mem.ind]
      co <- expand.grid(grid$xyCoords$y, grid$xyCoords$x)[2:1]
      le <- nrow(co)
      aux <- vapply(1:n.mem, FUN.VALUE = numeric(le), FUN = function(x) {
            z <- asub(grid[["Data"]], idx = x, dims = mem.ind, drop = TRUE)
            z <- unname(abind(z, along = -1L))
            attr(z, "dimensions") <- c("time", "lat", "lon")
            array3Dto2Dmat(z)
      })
      colnames(aux) <- paste0("Member_", 1:n.mem)
      aux <- data.frame(aux)
      co <- sp::coordinates(sp::points2grid(sp::SpatialPoints(co), tolerance = tolerance))
      co <- co[order(co[,1], co[,2]),]
      df <- cbind.data.frame(co, aux)
      sp::coordinates(df) <- c(1,2)
      sp::gridded(df) <- TRUE
      df <- as(df, "SpatialGridDataFrame")            
      # Argument list
      if (bt != "none") {
            uri <- switch(bt,
                          "coastline" = system.file("coastline.rda", package = "downscaleR"),
                          "countries" = system.file("countries.rda", package = "downscaleR"))
            load(uri)      
            if (is.null(arg.list[["sp.layout"]])) {
                  arg.list[["sp.layout"]] <- list(l1)
            } else {
                  arg.list[["sp.layout"]][[length(arg.list[["sp.layout"]]) + 1]] <- l1
            } 
      }
      if (is.null(arg.list[["col.regions"]])) {
            jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
            arg.list[["col.regions"]] <- jet.colors(31)
      }
      arg.list[["obj"]] <- df
      arg.list[["asp"]] <- 1
      print(do.call("spplot", arg.list))
}      


