#     filterGrid.R Apply time filters to a grid 
#
#     Copyright (C) 2016 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Time filtering
#' @description Apply a filter along the time dimension of a grid
#' @param grid Input grid (possibly multimember)
#' @template templateParallelParams
#' @param filter A vector of filter coefficients (see \code{filter for details})
#' @param method Either \code{"convolution"} or \code{"recursive"} (see \code{filter} for details)
#' @template templateParallel
#' @param sides For convolution filters only. If \code{sides = 1} the filter coefficients
#' are for past values only; if \code{sides = 2} they are centred around lag 0.
#' See \code{filter} for more details.
#' @param ... Further arguments passed to \code{filter}.
#' @return A time-filtered grid. 
#' @details  A wrapper of function \code{\link{filter}}
#' 
#' @export
#' @importFrom parallel stopCluster parApply
#' @importFrom stats filter
#' 
#' @author J Bedia
#'
#' @examples
#' data(iberia_ncep_ta850)
#' plot(iberia_ncep_ta850[["Data"]][,3,3], ty = 'l')
#' fgrid30 <- filterGrid(iberia_ncep_ta850, method = "convolution", filter = rep(1/30,30), sides = 1)
#' lines(fgrid30[["Data"]][,3,3], col = 'red')
#' fgrid90 <- filterGrid(iberia_ncep_ta850, method = "convolution", filter = rep(1/90,90), sides = 1)
#' lines(fgrid90[["Data"]][,3,3], col = 'green')
#' legend("top", c("raw","30 day MA", "90-day MA"), lty = 1, col = c(1,2,3), ncol = 3)

filterGrid <- function(grid, filter, method = c("convolution", "recursive"), sides = 1,
                       parallel = FALSE, max.ncores = 16, ncores = NULL, ...) {
      arg.list <- list(...)
      arg.list[["method"]] <- match.arg(method, choices = c("convolution", "recursive"))
      arg.list[["filter"]] <- filter
      arg.list[["sides"]] <- sides
      arr <- grid$Data
      refdim <- dim(arr)
      dimNames <- getDim(grid)
      mar <- grep("^time$", dimNames, invert = TRUE)
      stopifnot(!is.null(mar))
      ntimes <- dim(arr)[grep("^time$", dimNames)]
      stopifnot(ntimes > 1)
      message("[", Sys.time(), "] - Filtering ...")
      arr <- if (isTRUE(parallel)) {
                 parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
                 on.exit(parallel::stopCluster(parallel.pars$cl))
                 unname(parallel::parApply(cl = parallel.pars$cl, arr, MARGIN = mar, FUN = function(y) {
                     arg.list[["x"]] <- y
                     do.call("filter", arg.list)
                 }))
                 
             } else {
                 unname(apply(arr, MARGIN = mar, FUN = function(y) {
                     arg.list[["x"]] <- y
                     do.call("filter", arg.list)
                 }))
             }
      message("[", Sys.time(), "] - Done.")
      newdim <- dim(arr)
      if (!identical(newdim, refdim))  arr <- aperm(arr, perm = match(newdim, refdim))
      grid$Data <- arr
      attr(grid$Data, "dimensions") <- dimNames
      attr(grid$Variable, "filter:method") <- arg.list
      return(grid)
}
# End
