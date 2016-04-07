##     climatology.R Compute a grid climatology
##
##     Copyright (C) 2016 Santander Meteorology Group (http://www.meteo.unican.es)
##
##     This program is free software: you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 3 of the License, or
##     (at your option) any later version.
## 
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
## 
##     You should have received a copy of the GNU General Public License
##     along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' @title Compute a grid climatology
#' @description Calculates the climatology (i.e., complete temporal aggregation, typically the mean)
#' of the input grid.
#' @param grid Input grid
#' @param clim.fun Function to compute the climatology. This is specified as a list,
#'  indicating the name of the aggregation function in first place (as character), and other optional arguments
#'   to be passed to the aggregation function. Default to mean (i.e., \code{clim.fun = list(FUN="mean",na.rm = TRUE)}).
#' @param by.member Logical. In case of multimember grids, should the climatology be computed sepparately
#' for each member (\code{by.member=TRUE}), or a single climatology calculated from the ensemble mean
#'  (\code{by.member=FALSE})?. Default to \code{TRUE}.
#' @template templateParallelParams
#' @return A grid corresponding to the climatology. See details.
#' @template templateParallel
#' @details Two attributes are appended to the grid: 
#' \itemize{
#' \item \code{climatology:fun}, added to the \code{Data} component of the grid,
#' indicating the function used to compute the climatology.
#' \item \code{season}, added to the \code{Dates} component (if not yet existing), in order to provide information
#'  on the season for which the climatology has been computed.
#' }
#' @importFrom parallel stopCluster
#' @importFrom parallel parApply
#' @author J. Bedia
#' @export
#' @examples \dontrun{
#' #' # Maximum July surface temp forecast climatology
#' data("tasmax_forecast")
#' # (Parallelization option has no effect under WinOS)
#' # Aggregate all members before computing the climatology
#' tx_mean.clim <- climatology(tasmax_forecast,
#'                             by.member = FALSE,
#'                             parallel = TRUE)
#' # Note that time dimension is not dropped, and the new attributes
#' str(tx_mean.clim$Data)
#' str(tx_mean.clim$Dates)
#' # Compute a climatology for each member sepparately
#' tx_mean_15mem.clim <- climatology(tasmax_forecast,
#'                                   by.member = TRUE,
#'                                   parallel = TRUE)
#' str(tx_mean_15mem.clim$Data)
#' # 9 different climatologies, one for each member
#' 
#' # Flexible aggregation function definition:
#' # Example: climatology of the absolute maximum daily precipitation (Winter 2000):
#' data("NCEP_Iberia_tp")
#' tpmax <- climatology(NCEP_Iberia_tp,
#'                     clim.fun = list(FUN = "max"))
#' plotMeanGrid(tpmax)
#' }


# load("~/workspace/EUPORIAS/data/NCEP_psl_monthly_JA_1950_2010.rda", verbose = TRUE)
# load("~/workspace/EUPORIAS/data/S4_15_psl_monthly_JA_1981_2010_nnNCEPgrid.rda", verbose = TRUE)
# grid = psl.s4
# clim.fun = list(FUN = "mean", na.rm = TRUE)
# by.member = TRUE
# parallel = TRUE
# max.ncores = 16
# ncores = NULL

climatology <- function(grid,
                        clim.fun = list(FUN = "mean", na.rm = TRUE),
                        by.member = TRUE,
                        parallel = FALSE,
                        max.ncores = 16,
                        ncores = NULL) {
      parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
      dimNames <- attr(grid[["Data"]], "dimensions")
      ## Member aggregation 
      if ("member" %in% dimNames && !isTRUE(by.member)) {
            grid <- memberAggregation(grid,
                                      aggr.mem = list(FUN = "mean", na.rm = TRUE),
                                      parallel,
                                      max.ncores,
                                      ncores)
            dimNames <- attr(grid[["Data"]], "dimensions")
      }
      mar <- grep("^time$", dimNames, invert = TRUE)
      if (length(dimNames) == length(mar)) stop("Time dimension not found", call. = FALSE)
      arg.list <- c(clim.fun, list("MARGIN" = mar), list("X" = grid[["Data"]]))
      clim <- if (parallel.pars$hasparallel) {
            message("[", Sys.time(), "] - Computing climatology in parallel...")
            arg.list[["cl"]] <- parallel.pars$cl
            on.exit(parallel::stopCluster(parallel.pars$cl))
            do.call("parApply", arg.list)      
      } else {
            message("[", Sys.time(), "] - Computing climatology...")
            apply(grid[["Data"]], MARGIN = mar, FUN = "mean", na.rm = TRUE)      
            do.call("apply", arg.list)      
      }
      message("[", Sys.time(), "] - Done.")
      clim <- abind(clim, along = -1L)
      dimNames.aux <- c("time", dimNames[-grep("^time", dimNames)])
      attr(clim, "dimensions") <- dimNames.aux
      ## Dimension reordering
      perm <- na.omit(match(c("member","time","lat","lon"), dimNames.aux))
      clim <- aperm(clim,perm)
      grid[["Data"]] <- unname(clim)
      ## Attributes of Data
      attr(grid$Data, "dimensions") <- dimNames
      attr(grid$Data, "climatology:fun") <- clim.fun[["FUN"]]
      ## Date adjustment
      attr(grid$Dates, "season") <- getSeason(grid)
      grid$Dates$start <- grid$Dates$start[1]
      grid$Dates$end <- tail(grid$Dates$end, 1)
      return(grid)
}


