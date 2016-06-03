# detrendGrid.R Linear detrending of a grid
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

#' @title Linear detrending
#' @description Perform a linear detrending along the time dimension of a grid
#' @param grid Input grid (possibly multimember)
#' @template templateParallelParams
#' @return A detrended grid. 
#' @details  Performs a simple linear detrending by fitting a linear model and retaining the residuals.
#' An attribute indicating the linear detrending is added to the \code{Variable} component of the output grid.
#' 
#' In the presence of missing data in the time series, it operates by filtering them prior to linear model fitting. The
#' missing data positions are then restored back to the output detrended series.
#' 
#' @template templateParallel
#' @export
#' @importFrom parallel stopCluster parApply
#' @importFrom stats lm
#' @export
#' @author J Bedia
#' @examples 
#' data("iberia_ncep_ta850")
#' monthly <- aggregateGrid(iberia_ncep_ta850, aggr.m = list(FUN = "mean"))
#' plot(monthly$Data[,4,2], ty = 'l')
#' abline(reg = lm(monthly$Data[,4,2] ~ I(1:length(monthly$Data[,4,2]))))
#' det <- detrendGrid(monthly, parallel = FALSE)
#' # Detrended series in red
#' lines(det$Data[,4,2], col = "red")
#' abline(reg = lm(det$Data[,4,2] ~ I(1:length(det$Data[,4,2]))), col = "red")

detrendGrid <- function(grid, parallel = FALSE, max.ncores = 16, ncores = NULL) {
      arr <- grid$Data
      refdim <- dim(arr)
      dimNames <- getDim(grid)
      mar <- grep("^time$", dimNames, invert = TRUE)
      ntimes <- dim(arr)[grep("^time$", dimNames)]
      x <- 1:ntimes
      message("[", Sys.time(), "] - Detrending...")
      arr <- if (isTRUE(parallel)) {
            parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
            on.exit(parallel::stopCluster(parallel.pars$cl))
                  unname(parallel::parApply(cl = parallel.pars$cl, arr, MARGIN = mar, FUN = function(y) {
                  out <- rep(NA, ntimes)
                  ind <- intersect(which(!is.na(y)), x)
                  out[ind] <- tryCatch(expr = summary(lm(y ~ I(x)))$resid + mean(y, na.rm = TRUE),
                                       error = function(err) {
                                             out
                                       })
                  return(out)
                  })
                  
            )
      } else {
            unname(apply(arr, MARGIN = mar, FUN = function(y) {
                  out <- rep(NA, ntimes)
                  ind <- intersect(which(!is.na(y)), x)
                  out[ind] <- tryCatch(expr = summary(lm(y ~ I(x), subset = ind))$resid + mean(y, na.rm = TRUE),
                                       error = function(err) {
                                             out
                                       })
                  return(out)
                  })
            )
      }
      message("[", Sys.time(), "] - Done.")
      newdim <- dim(arr)
      if (!identical(newdim, refdim))  arr <- aperm(arr, perm = match(newdim, refdim))
      grid$Data <- arr
      attr(grid$Data, "dimensions") <- dimNames
      attr(grid$Variable, "detrended:method") <- "linear"
      return(grid)
}
# End
