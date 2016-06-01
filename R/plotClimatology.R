#' @title Lattice plot methods for climatological grids
#' @description A wrapper for the lattice (trellis) plot methods for spatial data in \code{sp::spplot}
#' @param grid Input grid
#' @param tolerance Precision, used to which extent points are exactly on a grid. See details.
#' @param backdrop.theme Reference geographical lines to be added to the plot. See Details. 
#' @param ... Further arguments passed to \code{spplot}
#' @details The function applies the \code{\link[sp]{spplot}} method after conversion of the climatological map(s) to a
#'  \code{SpatialGridDataFrame}.
#'  
#'  \strong{Multigrids}
#'  
#'  Multigrids of climatologies can be created using \code{makeMultiGrid} 
#'  for trellis visualization of different variables, or for instance, for the comparison of
#'  raw and corrected/downscaled scenarios side to side. In case of multimember multigrids, 
#'  the function will internally compute the ensemble mean of each variable in the multigrid
#'   for representation (with a message).
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
#' 
#' @export
#' 
#' @author J. Bedia
#' @seealso \code{\link{climatology}}. See \code{\link[sp]{spplot}} in package \pkg{sp} for further information on
#' plotting capabilities and options
#' @examples \donttest{
#' data(tasmax_forecast)
#' # Climatology is computed:
#' clim <- climatology(tasmax_forecast, by.member = TRUE)
#' # Some grids not perfectly regular (e.g. CFSv2) may yield a spatial tolerance error:
#' try(plotClimatology(clim))
#' # A tolerance value is suggested:
#' plotClimatology(clim, tolerance = 0.00056)
#' # Geographical lines can be added using the argument 'backdrop.theme':
#' plotClimatology(clim, tolerance = 0.00056, backdrop.theme = "coastline")
#' 
#' # Further arguments can be passed to 'spplot'...
#' 
#' # ... a subset of members to be displayed, using 'zcol':
#' plotClimatology(clim,
#'                 tolerance = 0.00056,
#'                 backdrop.theme = "coastline",
#'                 zcol = 1:4)
#'                 
#' # ... regional focuses (e.g. the Iberian Peninsula):
#' plotClimatology(clim,
#'                 tolerance = 0.00056,
#'                 backdrop.theme = "countries",
#'                 xlim = c(-10,5), ylim = c(35,44),
#'                 zcol = 1:4,
#'                 scales = list(draw = TRUE))
#' 
#' # Changing the default color palette and ranges:
#' plotClimatology(clim,
#'                 tolerance = 0.00056,
#'                 backdrop.theme = "coastline",
#'                 zcol = 1:4,
#'                 col.regions = heat.colors(27), at = seq(10,37,1))
#'                 
#' # For ensemble means climatology should be called with 'by.member' set to FALSE:
#' clim <- climatology(tasmax_forecast, by.member = FALSE)
#' 
#' # Adding contours to the plot is direct with argument 'contour':
#' 
#' plotClimatology(clim,
#'                 tolerance = 0.00056,
#'                 scales = list(draw = TRUE),
#'                 contour = TRUE,
#'                 main = "tasmax Predictions July Ensemble Mean")
#' }

plotClimatology <- function(grid,
                            tolerance = sqrt(.Machine$double.eps),
                            backdrop.theme = "none",
                            ...) {
      arg.list <- list(...)
      bt <- match.arg(backdrop.theme, choices = c("none", "coastline", "countries"))
      dimNames <- getDim(grid)
      ## Multigrids are treated as realizations, previously aggregated by members if present
      is.multigrid <- "var" %in% dimNames
      if (is.multigrid) {
            if ("member" %in% dimNames) {
                  mem.ind <- grep("member", dimNames)
                  n.mem <- dim(grid[["Data"]])[mem.ind]
                  if (n.mem > 1) message("NOTE: The multimember mean will be displayed for each variable in the multigrid")
                  grid <- aggregateGrid(grid, aggr.mem = list(FUN = "mean", na.rm = TRUE))
                  dimNames <- getDim(grid)
            }
            attr(grid[["Data"]], "dimensions") <- gsub("var", "member", dimNames)      
      }
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
      if (is.multigrid) {
            vname <- attr(grid$Variable, "longname")
            if (!is.null(grid$Variable$level)) {
                  auxstr <- paste(vname, grid$Variable$level, sep = "@")
                  vname <- gsub("@NA", "", auxstr)
            }
            vname <- gsub("\\s", "_", vname)
      } else {
            vname <- paste0("Member_", 1:n.mem)
      }
      colnames(aux) <- vname
      aux <- data.frame(aux)
      co <- sp::coordinates(sp::points2grid(sp::SpatialPoints(co), tolerance = tolerance))
      co <- co[order(co[,1], co[,2]),]
      df <- cbind.data.frame(co, aux)
      sp::coordinates(df) <- c(1,2)
      sp::gridded(df) <- TRUE
      df <- as(df, "SpatialGridDataFrame")            
      ## Backdrop theme -----
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
      ## Default colorbar -------
      if (is.null(arg.list[["col.regions"]])) {
            jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
            arg.list[["col.regions"]] <- jet.colors(31)
      }
      ## Other args --------
      arg.list[["obj"]] <- df
      arg.list[["asp"]] <- 1
      print(do.call("spplot", arg.list))
}      


