# downscale.R Perfect-prog downscaling methods
#
#     Copyright (C) 2015 Santander Meteorology Group (http://www.meteo.unican.es)
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
#' @param n.neigh Applies only when \code{method="analogs"} (otherwise ignored). Integer indicating the number of closest neigbours to retain for analog construction. Default to 1.
#' @param sel.fun Applies only when \code{method="analogs"} (otherwise ignored). Criterion for the construction of analogs when several neigbours are chosen. Ignored when \code{n.neig = 1}.
#' Current values are \code{"random"} (the default) and \code{"mean"}. See details.
#' @param analog.dates Logical flag indicating whether the dates of the analogs should be returned. See the analogs section.
#' @param pr.threshold Value below which precipitation amount is considered zero 
#' @param n.pcs Integer indicating the number of EOFs to be used as predictors
#' @param cross.val Should cross-validation be performed?
#' @template templateParallelParams 
#' 
#' @details
#' 
#' \strong{Spatial consistency}
#' 
#' Several checks of spatial consistency are performed. In particular, note that both \code{pred} (reanalysis) and \code{sim}
#'  (model simulations) should be in the same grid. This consistency must be ensured by the user prior to entering these arguments,
#' for instance by means of the \code{\link{interpData}} function in conjunction with the \code{\link{getGrid}} method.
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
#' \code{n.neigh} > 1). In this case the \code{"random"} choice randomly selects one of the \code{n.neigh} neighbours,
#'  while the \code{"mean"} choice will compute their average.
#' 
#' \strong{Analog dates}
#' 
#' If the argument \code{analog.dates} is set to TRUE, these will be returned as a global attribute named \code{"analog.dates"}.
#' Note that the analog dates can be only returned for one single neighbour selections (i.e. \code{n.neigh = 1}),
#'  otherwise giving an error. The analog dates are different for each member in case of 
#'  multimember downscaling, and are returned as a list, each element of the list corresponding to one member. 
#' 
#' @template templateParallel
#' 
#'   
#' @seealso \code{\link{prinComp}} for details on principal component/EOF analysis,
#' \code{rescaleMonthlyMeans} for data pre-processing,
#' \code{\link{makeMultiField}} for multifield creation
#' \code{\link{loadGridData}} and \code{\link{loadStationData}} for loading fields and station data respectively.
#' 
#' @export 
#' @family downscaling
#' @author J Bedia and M Iturbide

# load("ignore/juaco/data/obsPredSim_NCEP.Rdata", verbose = TRUE)
# 
# plotMeanField(pred)
# 
# 
# # # 
# pca.pred <- prinComp(pred, v.exp = .99)
# str(pca.pred)

# str(sim)

# str(a)
# 
# str(pred)
# str(pca.pred)
# # 
# # a <- ppModelSetup(obs.precip, pca.pred, sim)
# # str(a)
# # str(pca.pred)
# str(pred)
# str(modelPars)
# parallel = TRUE
# ncores = NULL
# max.ncores = 16
# 
# 
# a <- downscale(obs.tmean, pred, sim, method = "analogs", parallel = TRUE)
# 
# b <- downscale(obs.tmean, pred, sim, method = "analogs", parallel = FALSE)


downscale <- function(obs,
                      pred,
                      sim = NULL,
                      method = c("analogs", "glm"),
                      n.neigh = 1,
                      sel.fun = c("random", "mean"),
                      analog.dates = FALSE,
                      pr.threshold = .1,
                      n.pcs = NULL,
                      cross.val = c("none", "loocv"),
                      parallel = FALSE,
                      max.ncores = 16,
                      ncores = NULL) {
      cross.val <- match.arg(cross.val, choices = c("none", "loocv"))
      parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
      modelPars <- ppModelSetup(obs, pred, sim)
      if (cross.val == "none") {
            switch(method,
                  "analogs" = analogs(obs = obs,
                                      modelPars = modelPars,
                                      n.neigh = n.neigh,
                                      sel.fun = sel.fun,
                                      analog.dates = analog.dates,
                                      parallel.pars = parallel.pars),
                  "glm" = glimpr(obs = obs,
                                 modelPars = modelPars,
                                 pr.threshold = pr.threshold,
                                 n.pcs = n.pcs)
            )
      } else if (cross.val == "loocv") {
            if (!is.null(sim) & !identical(pred, sim)) {
                  
                  
            }
      }
}
# identical(pred,sim)
# # ppModelSetup  la matriz pred tiene que tener un atributo que indica que es un producto de PCA o no
# a <- downscale(obs = obs.precip, pred = pred, sim = sim, method = "analogs")
# str(a)
# # 
# # range(obs.precip$Data)
# obs = obs.precip
# str(modelPars)
