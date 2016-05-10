##     rescaleGrid.R Grid rescaling
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


#' @title Grid rescaling
#' @description Rescale a grid
#' @param grid Input grid to be rescaled
#' @param grid.clim Reference climatology to be subtracted to the input grid. If \code{NULL} (the default),
#' the climatology is directly computed from the grid via \code{\link{climatology}}, either member-by-member
#' or using the ensemble mean climatology, as specified by the \code{by.member} flag.
#' @param ref Reference grid. After subtracting to grid its climatology (as defined by \code{grid.clim}), the mean
#' climatology of this additional grid is added. Default to NULL, so no reference grid is used.
#' @param by.member Logical. In case of multimember grids, should the climatology be computed sepparately
#' for each member (\code{by.member=TRUE}), or a single climatology calculated from the ensemble mean
#'  (\code{by.member=FALSE})?. Default to \code{TRUE}. Argument passed to \code{\link{climatology}}.
#' @template templateParallelParams
#' @details The reference grid (\code{ref}) is used to correct the input grid, as follows:
#' 
#' \deqn{grid' = grid - mu_grid + mu_ref}
#' 
#' , where \emph{mu} corresponds to the monthly climatological mean considering the training period,
#' and \emph{grid'} is the corrected test data. The way \emph{mu_ref} and \emph{mu_grid} are computed in case
#' of multimember grids is controlled by the argument \code{by.member}.
#' 
#' The \code{ref} usually corresponds to the control run of the GCM in the training period in climate change applications,
#' or the hindcast data for the training period in s2d applications. Note that by default \code{ref = NULL}. In this 
#' case it will be assumed to be the \code{pred} grid. This can be used for instance when train and test correspond
#' to the same model.
#' @importFrom abind abind asub
#' @importFrom stats na.omit
#' @return A rescaled grid
#' @export



# load("~/workspace/EUPORIAS/data/NCEP_psl_monthly_JA_1950_2010.rda", verbose = TRUE)
# load("~/workspace/EUPORIAS/data/S4_15_psl_monthly_JA_1981_2010_nnNCEPgrid.rda", verbose = TRUE)
# plotMeanGrid(psl)
# plotMeanGrid(psl.s4)
# psl.s4.resc <- rescaleGrid(psl.s4, grid.clim = NULL, ref = psl, parallel =F, by.member = FALSE)
# 
# plotMeanGrid(psl.s4.resc)
# 
# 
# ref <- subsetGrid(psl, years = 1981:2010)
# grid <- psl.s4
# grid.clim = NULL
# by.member = FALSE
# parallel = FALSE
# max.ncores = 16
# ncores = NULL



rescaleGrid <- function(grid,
                        grid.clim = NULL,
                        ref = NULL,
                        by.member = TRUE,
                        parallel = FALSE,
                        max.ncores = 16,
                        ncores = NULL) {
      parallel.pars <- parallelCheck(parallel, max.ncores, ncores) 
      if (is.null(grid.clim)) {
            grid.clim <- climatology(grid,
                                     by.member = by.member,
                                     parallel = parallel,
                                     max.ncores = max.ncores,
                                     ncores = ncores)
      }
      if (!identical(getSeason(grid.clim), getSeason(grid))) {
            stop("Seasons of input grid and grid.clim do not match")
      }
      if (!is.null(ref) && !identical(getSeason(grid), getSeason(ref))) {
            stop("Seasons of input grid and reference grid do not match")
      }
      if (!is.null(ref)) {
            ref.clim <- climatology(ref,
                                    by.member = by.member,
                                    parallel = parallel,
                                    max.ncores = max.ncores,
                                    ncores = ncores)
            
            if ("member" %in% getDim(grid.clim) && !("member" %in% getDim(ref))) {
                  n.mem <- dim(grid.clim[["Data"]])[grep("member", getDim(grid.clim))]
                  aux <- rep(list(ref.clim[["Data"]]), n.mem)
                  dimNames <- getDim(ref.clim)
                  ref.clim[["Data"]] <- unname(do.call("abind", c(aux, along = -1L)))
                  aux <- NULL
                  attr(ref.clim[["Data"]], "dimensions") <- c("member", dimNames)
            }
      } else {
            ref.clim <- list()
            ref.clim[["Data"]] <- array(data = 0, dim = dim(grid.clim[["Data"]]))
            attr(ref.clim[["Data"]], "dimensions") <- attr(grid.clim[["Data"]], "dimensions")
      }
      clim <- grid[["Data"]]
      dimNames <- getDim(grid)
      ind.time <- grep("^time", dimNames)
      n.times <- dim(clim)[ind.time]
      Xc <- drop(grid.clim[["Data"]])
      Xref <- drop(ref.clim[["Data"]])
      aux.list <- lapply(1:n.times, function(x) {
            X <- asub(clim, idx = x, dims = ind.time)
            X - Xc + Xref
      })
      clim <- do.call("abind", c(aux.list, along = -1L))
      dimNames.aux <- c("time", dimNames[-grep("^time", dimNames)])
      attr(clim, "dimensions") <- dimNames.aux
      ## Dimension reordering
      perm <- na.omit(match(c("member","time","lat","lon"), dimNames.aux))
      clim <- aperm(clim, perm)
      grid[["Data"]] <- unname(clim)
      attr(grid[["Data"]], "dimensions") <- dimNames
      return(grid)
}

