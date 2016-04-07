# aggregateGrid.R Grid aggregation along selected dimensions
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


#' @title Grid aggregation along selected dimensions
#' @description Aggregates a grid along the target dimensions through aggregation function specification.
#' @param grid a grid or multigrid to be aggregated.
#' @param aggr.d Daily aggregation function (for sub-daily data only). A list indicating the name of the
#'  aggregation function in first place, and other optional arguments to be passed to the aggregation function. See the examples.
#' @param aggr.m Same as \code{aggr.d}, but indicating the monthly aggregation function. 
#' @param aggr.y Same as \code{aggr.d}, but indicating the annual aggregation function. 
#' @param aggr.mem Same as \code{aggr.d}, but indicating the function for computinh the member aggregation.
#' @param aggr.lat Same as \code{aggr.d}, indicating the aggregation function to be applied along latitude.
#' @param aggr.lon Same as \code{aggr.lat}, but for longitude.
#' @template templateParallelParams
#' @return A grid or multigrid aggregated along the chosen dimension(s).
#' @details
#' 
#' \strong{Aggregation function definition}
#' 
#' The aggregation functions are specified in the form of a named list of the type \code{FUN = "function", ...}, where
#' \code{...} are further arguments passes to FUN. This allows for a flexible definition of aggregation functions, that are 
#' internally passes to \code{\link{tapply}}. Note that the name of the function is indicated as a character string.
#' 
#' \strong{Member aggregation}
#' 
#' The function preserves the metadadata associated with member information (i.e. initialization dates and member names) after
#' aggregation. In addition, an attribute indicating the member aggregation function is added to the \code{Variable} component.
#' 
#' 
#' \strong{Temporal aggregation}
#'  
#' To annually or monthly aggregate data, \code{aggr.d} and/or \code{aggr.m} functions are specified.
#' Aggregatikons need to be specified from bottom to top, so for instance, if the data in the grid is sub-daily
#' and \code{aggr.d} is not specified, an error will be given for monthly or annual aggregation requests. Similarly,
#' annual aggregations require a previous specification of daily and monthly aggregation, when applicable. Special attributes
#' in the \code{Variable} component indicate the aggregation undertaken.
#' 
#' @template templateParallel
#' @author M. Iturbide, M. de Felice, J. Bedia 
#' @export
#' @examples \dontrun{
#' data("iberia_tasmax")
#' ## Aggregating members
#' # Ensemble mean
#' mn <- aggregateGrid(grid = tasmax_forecast, aggr.mem = list("mean", na.rm = TRUE))
#' # Ensemble 90th percentile
#' ens90 <- aggregateGrid(grid = tasmax_forecast,
#'                        aggr.mem = list("quantile", probs = 0.9, na.rm = TRUE))
#' par(mfrow = c(1,2))
#' plotMeanGrid(mn)
#' plotMeanGrid(ens90)
#' par(mfrow = c(1,1))
#' 
#' ## Monthly aggregation
#' monthly.mean <- aggregateGrid(tasmax_forecast, aggr.m = list(FUN = mean, na.rm = TRUE))
#' 
#' ## Several dimensions ca be aggregated in one go:
#' mm.mean <- aggregateGrid(tasmax_forecast,
#'                aggr.mem = list(FUN = "mean", na.rm = TRUE),
#'                aggr.m = list(FUN = "mean", na.rm = TRUE))
#' }

aggregateGrid <- function(grid,
                          aggr.mem = list(FUN = NULL),
                          aggr.d = list(FUN = NULL),
                          aggr.m = list(FUN = NULL),
                          aggr.y = list(FUN = NULL),
                          aggr.lat = list(FUN = NULL),
                          aggr.lon = aggr.lat,
                          parallel = FALSE,
                          max.ncores = 16,
                          ncores = NULL) {
      if (!is.null(aggr.mem$FUN)) {
            grid <- memberAggregation(grid, aggr.mem, parallel, max.ncores, ncores)
      }
      if (!is.null(aggr.d$FUN)) {
            grid <- timeAggregation(grid, "DD", aggr.d, parallel, max.ncores, ncores)
      }
      if (!is.null(aggr.m$FUN)) {
            grid <- timeAggregation(grid, "MM", aggr.m, parallel, max.ncores, ncores)
      }
      if (!is.null(aggr.y$FUN)) {
            grid <- timeAggregation(grid, "YY", aggr.y, parallel, max.ncores, ncores)
      }
      return(grid)
}

#' @title Member aggregation
#' @description Aggregate a grid along its member dimension
#' @param grid A multimember grid to apply the aggregation
#' @param aggr.mem Character string indicatins the aggregation function
#' @param parallel.pars Arguments defining the parallelization options, as passed by \code{\link{parallelCheck}}
#' @details The function preserves the metadadata associated with member information (i.e. initialization dates and member names). In addition,
#' an attribute indicating the member aggregation function is added to the \code{Variable} component.
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster
#' @keywords internal
#' @author J. Bedia

memberAggregation <- function(grid, aggr.mem, parallel, max.ncores, ncores) {
      dimNames <- attr(grid$Data, "dimensions")
      if (!"member" %in% dimNames) {
            message("Not a multimember grid: 'aggr.mem' option was ignored.")
      } else {
            parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
            attr.all <- attributes(grid$Data)
            mar <- grep("member", dimNames, invert = TRUE)
            attr.all$dim <- attr.all$dim[mar]
            attr.all$dimensions <- dimNames[mar]
            aggr.mem[["MARGIN"]] <- mar
            aggr.mem[["X"]] <- grid$Data
            out <- if (parallel.pars$hasparallel) {
                  message("[", Sys.time(), "] - Aggregating members in parallel...")
                  on.exit(parallel::stopCluster(parallel.pars$cl))
                  aggr.mem[["cl"]] <- parallel.pars$cl
                  do.call("parApply", aggr.mem)
            } else {
                  message("[", Sys.time(), "] - Aggregating members...")
                  do.call("apply", aggr.mem)
            }
            message("[", Sys.time(), "] - Done.")
            grid[["Data"]] <- out
            if (any(names(attr.all) != "dim" & names(attr.all) != "dimensions")) {
                  attributes(grid$Data) <- attr.all[grep("^dim$|^dimensions$", names(attr.all), invert = TRUE)]
            }
            dimNames <- dimNames[-grep("member", dimNames)]
            attr(grid$Data, "dimensions") <- dimNames
            attr(grid$Variable, "member_agg_cellfun") <- aggr.mem[[1]]
      }
      return(grid)
}


#' @title Time aggregation
#' @description Aggregate a grid along its time dimension
#' @param grid A multimember grid to apply the aggregation
#' @param aggr.type Character string indicating the type of temporal aggregation: 
#' daily (\code{"DD"}), monthly (\code{"MM"}) or annual (\code{"YY"}).
#' @param aggr.fun One of \code{aggr.d}, \code{aggr.m} or \code{aggr.y} arguments, as passed by \code{aggregateGrid}
#' @param parallel.pars Arguments defining the parallelization options, as passed by \code{\link{parallelCheck}}
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster
#' @keywords internal
#' @author J. Bedia, M. Iturbide, M. de Felice


timeAggregation <- function(grid, aggr.type = c("DD","MM","YY"), aggr.fun, parallel, max.ncores, ncores) {
      aux.dates <- if ("var" %in% attr(grid$Data, "dimensions")) {
            grid$Dates[[1]]$start
      } else {
            grid$Dates$start
      }
      dff <- abs(difftime(aux.dates[1], aux.dates[2], units = "hours"))
      if (aggr.type == "DD" & dff >= 24) {
            message("Data is already daily: 'aggr.d' option was ignored.")
      } else if (aggr.type == "MM" & dff >= 672) {
            message("Data is already monthly: 'aggr.m' option was ignored.")
      } else if (aggr.type == "YY" & dff >= 8640) {
            message("Data is already 6annual: 'aggr.y' option was ignored.")
      } else {
            dimNames <- attr(grid$Data, "dimensions")
            # Attributes
            attr.all <- attributes(grid$Data)
            mar <- grep("^time", dimNames, invert = TRUE)
            day <- substr(aux.dates,9,10)
            mon <- substr(aux.dates,6,7)
            yr <- getYearsAsINDEX(grid)
            fac <- switch(aggr.type,
                          "DD" = paste0(yr,mon,day),
                          "MM" = paste0(yr,mon),
                          "YY" = yr)
            day <- mon <- yr <- aux.dates <- NULL
            arg.list <- c(aggr.fun, list("INDEX" = fac))
            type <- switch(aggr.type,
                          "DD" = "daily",
                          "MM" = "monthly",
                          "YY" = "annual")
            parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
            arr <- if (parallel.pars$hasparallel) {
                  message("[", Sys.time(), "] Performing ", type, " aggregation in parallel...")
                  on.exit(parallel::stopCluster(parallel.pars$cl))
                  parallel::parApply(cl = parallel.pars$cl, grid$Data, MARGIN = mar, FUN = function(x) {
                        arg.list[["X"]] <- x
                        do.call("tapply", arg.list)
                  })
            } else {
                  message("[", Sys.time(), "] Performing ", type, " aggregation...")
                  apply(grid$Data, MARGIN = mar, FUN = function(x) {
                        arg.list[["X"]] <- x
                        do.call("tapply", arg.list)
                  })
            }
            message("[", Sys.time(), "] Done.")
            # Array attributes -----------------
            if (length(dim(arr)) != length(dimNames)) arr <- abind(arr, along = -1) # Preserve time dimension if lost
            if (grep("^time", dimNames) > 1) {
                  arr <- aperm(arr, c(grep("^time", dimNames), grep("^time", dimNames, invert = TRUE)))     
            }
            grid$Data <- unname(arr)
            if (any(names(attr.all) != "dim" & names(attr.all) != "dimensions")) {
                  attributes(grid$Data) <- attr.all[grep("^dim$|^dimensions$", names(attr.all), invert = TRUE)]
            }
            attr(grid$Data, "dimensions") <- dimNames
            # Date adjustment ------------------
            if ("var" %in% dimNames) {
                  grid$Dates <- lapply(1:length(grid$Dates), function(x) {
                        list("start" = unname(tapply(grid$Dates[[x]]$start, INDEX = fac, FUN = min)),
                             "end" = unname(tapply(grid$Dates[[x]]$end, INDEX = fac, FUN = max)))
                  })
            } else {
                  grid$Dates <- list("start" = unname(tapply(grid$Dates$start, INDEX = fac, FUN = min)),
                                     "end" = unname(tapply(grid$Dates$end, INDEX = fac, FUN = max)))
            }
            # Temporal aggregation attributes --------
            attr(grid$Variable, paste0(type,"_agg_cellfun")) <- arg.list$FUN
      }
      return(grid)
}
                  
                  
                  

