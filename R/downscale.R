# downscale.R Perfect-prog downscaling methods
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

#' @title Perfect-prog downscaling
#' 
#' @description Workhorse function to call the different perfect-prog downscaling methods
#' 
#' @template templateObsPredSim
#' @param method Downscaling method
#' @param simulate Logical. Default is TRUE.
#' @param n.analogs Applies only when \code{method="analogs"} (otherwise ignored). Integer indicating the number of closest neigbours to retain for analog construction. Default to 1.
#' @param sel.fun Applies only when \code{method="analogs"} (otherwise ignored). Criterion for the construction of analogs when several neigbours are chosen. Ignored when \code{n.neig = 1}.
#' Current values are \code{"random"} (the default) and \code{"mean"}. See details.
#' @param analog.dates Logical flag indicating whether the dates of the analogs should be returned. See the analogs section.
#' @param wet.threshold Value below which precipitation amount is considered zero 
#' @param n.pcs Integer indicating the number of EOFs to be used as predictors
#' @param cross.val Should cross-validation be performed? methods available are leave-one-out ("loocv") and k-fold ("kfold"). The default 
#' option ("none") does not perform cross-validation.
#' @param folds Only requiered if cross.val = "kfold". A list of vectors, each containing the years to be grouped in 
#' the corresponding fold.
#' @template templateParallelParams 
#' 
#' @details
#' 
#' \strong{Spatial consistency}
#' 
#' Several checks of spatial consistency are performed. In particular, note that both \code{x} (reanalysis) and \code{newdata}
#'  (model simulations) should be in the same grid. This consistency must be ensured by the user prior to entering these arguments,
#' for instance by means of the \code{\link{interpGrid}} function in conjunction with the \code{\link{getGrid}} method.
#' 
#' \strong{Scaling and centering}
#' 
#' When the climate variables are used as predictors instead of the PCs, these are previously centered and scaled
#' using the mean and sigma parameters globally computed for the whole spatial domain (This is equivalent to the \dQuote{field})
#' method in the \code{\link{prinComp}} function. The simulation data will use the parameters obtained when scaling and centering
#' the predictors dataset. In case that the predictors come from a PC analysis object (as returned by \code{\link{prinComp}}), the
#' parameters for rescaling the simulation data are passed by the predictors.
#' 
#' 
#' @section Analogs:
#' 
#' \strong{Construction of analogs using multiple neighbours}
#' 
#' The argument \code{sel.fun} controls how the analogs are constructed when considering more than the first neighbour (argument
#' \code{n.analogs} > 1). In this case the \code{"random"} choice randomly selects one of the \code{n.analogs} neighbours,
#'  while the \code{"mean"} choice will compute their average.
#' 
#' \strong{Analog dates}
#' 
#' If the argument \code{analog.dates} is set to TRUE, these will be returned as a global attribute named \code{"analog.dates"}.
#' Note that the analog dates can be only returned for one single neighbour analog selections (i.e. \code{n.analogs = 1}),
#'  otherwise giving an error. The analog dates are different for each member in case of 
#'  multimember downscaling, and are returned as a list, each element of the list corresponding to one member. 
#' @template templateParallel
#' @seealso \code{\link{prinComp}} for details on principal component/EOF analysis,
#' \code{rescaleMonthlyMeans} for data pre-processing,
#' \code{\link{makeMultiGrid}} for multigrid construction
#' \code{loadGridData} and \code{loadStationData}, from package \pkg{loadeR}, for loading grids and station data respectively.
#' @export 
#' @family downscaling
#' @author J Bedia and M Iturbide

downscale <- function(y,
                      x,
                      newdata = NULL,
                      method = c("analogs", "glm"),
                      simulate = TRUE,
                      n.analogs = 1,
                      sel.fun = c("random", "mean"),
                      analog.dates = FALSE,
                      wet.threshold = .1,
                      n.pcs = NULL,
                      cross.val = c("none", "loocv", "kfold"),
                      folds = NULL,
                      parallel = FALSE,
                      max.ncores = 16,
                      ncores = NULL) {
      cross.val <- match.arg(cross.val, choices = c("none", "loocv", "kfold"))
      parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
      modelPars <- ppModelSetup(y, x, newdata)
      obs.orig <- y
      if ("station" %in% attr(y$Data, "dimensions")){
            stations <-  TRUE
      } else {
            stations <- FALSE
      }
      if (cross.val == "none") {
            down <- switch(method,
                           "analogs" = analogs(y = y,
                                               modelPars = modelPars,
                                               n.analogs = n.analogs,
                                               sel.fun = sel.fun,
                                               analog.dates = analog.dates,
                                               parallel.pars = parallel.pars),
                           "glm" = glimpr(y = y,
                                          modelPars = modelPars,
                                          wet.threshold = wet.threshold,
                                          n.pcs = n.pcs,
                                          simulate = simulate)
            )
      } else {
            if (!is.null(newdata)) {
                  message("'newdata' will be ignored for cross-validation")
            }
            if ("scaled:method" %in% names(attributes(x))) {
                  x$Dates$start <- attr(x, "dates_start")
            }
            years <- getYearsAsINDEX(x)
            modelPars.orig <- modelPars
            message("[", Sys.time(), "] Fitting models...")
            if (cross.val == "loocv") {
                  downi <- lapply(1:length(unique(years)), function(i) {
                        year.ind <- which(years == unique(years)[i])
                        modelPars$sim.mat[[1]] <- modelPars.orig$pred.mat[year.ind,]
                        modelPars$pred.mat <- modelPars.orig$pred.mat[-year.ind,]
                        if (stations == TRUE) {
                              y$Data <- obs.orig$Data[-year.ind,]
                        } else {
                              y$Data <- obs.orig$Data[-year.ind,,]
                        }
                        attr(y$Data, "dimensions") <- attr(obs.orig$Data, "dimensions")
                        message("Validation ", i, ", ", length(unique(years)) - i, " remaining")
                              if (method == "analogs") {
                                    suppressMessages(
                                    analogs(y = y,
                                    modelPars = modelPars,
                                    n.analogs = n.analogs,
                                    sel.fun = sel.fun,
                                    analog.dates = analog.dates,
                                    parallel.pars = parallel.pars)
                                    )
                              } else if (method == "glm") {
                                    suppressMessages(
                                    glimpr(y = y,
                                    modelPars = modelPars,
                                    wet.threshold = wet.threshold,
                                    n.pcs = n.pcs,
                                    simulate = simulate)
                                    )
                              }
                        }
                  )
                  down <- unname(do.call("abind", c(downi, along = 1)))
            } else if (cross.val == "kfold") {
                  if (is.null(folds)) {
                        stop("Fold specification is missing, with no default", call. = FALSE)
                  }
                  downi <- lapply(1:length(folds), function(i) {
                        year.ind <- lapply(1:length(folds[[i]]), function(k){
                              indstep <- which(years == folds[[i]][k])
                              return(indstep)
                        })
                        year.ind <- unname(abind(year.ind, along = 1))
                        modelPars$sim.mat[[1]] <- modelPars.orig$pred.mat[year.ind,]
                        modelPars$pred.mat <- modelPars.orig$pred.mat[-year.ind,]
                        if (stations == TRUE) {
                              y$Data <- obs.orig$Data[-year.ind,]
                        } else {
                              y$Data <- obs.orig$Data[-year.ind,,]
                        }
                        attr(y$Data, "dimensions") <- attr(obs.orig$Data, "dimensions")
                        message("Validation ", i, ", ", length(folds) - i, " remaining")
                              if (method == "analogs") {
                                    suppressMessages(
                                          analogs(y = y,
                                                modelPars = modelPars,
                                                n.analogs = n.analogs,
                                                sel.fun = sel.fun,
                                                analog.dates = analog.dates,
                                                parallel.pars = parallel.pars)
                                    )
                              } else if (method == "glm") {
                                    suppressMessages(
                                          glimpr(y = y,
                                                modelPars = modelPars,
                                                wet.threshold = wet.threshold,
                                                n.pcs = n.pcs)
                                    )
                              }
                        })
                  down <- unname(do.call("abind", c(downi, along = 1)))
            }
      }
      # Data array - rename dims
      dimNames <- renameDims(obs.orig, modelPars$multi.member)
      obs.orig$Data <- down
      attr(obs.orig$Data, "dimensions") <- dimNames
      attr(obs.orig$Data, "downscaling:method") <- method
      attr(obs.orig$Data, "downscaling:cross-validation") <- cross.val
      attr(obs.orig$Data, "downscaling:simulation_data") <- modelPars$sim.dataset
      attr(obs.orig$Data, "downscaling:n_pcs") <- n.pcs
      # Date replacement
      obs.orig$Dates <- modelPars$sim.dates 
      message("[", Sys.time(), "] Done.")
      return(obs.orig)
}

