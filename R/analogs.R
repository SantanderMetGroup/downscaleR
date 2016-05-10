#' @title Analog downscaling
#' @description Implementation of the downscaling analogs method
#' @param y A grid or station data containing the observed climate data for the training period
#' @param modelPars Output object from function ppModelSetup containing the predictors and test data
#' @param n.analogs Integer indicating the number of closest neigbours to retain for analog construction. Default to 1.
#' @param sel.fun Criterion for the construction of analogs when several neigbours are chosen. Ignored when \code{n.neig = 1}.
#' Current values are \code{"random"} (the default) and \code{"mean"}. See details.
#' @param analog.dates Logical flag indicating whether the dates of the analogs should be returned. If set to TRUE,
#'  the analog dates will be returned as a global attribute named \code{"analog.dates"}. The analog
#'   dates can be only returned for one single neighbour selections (argument \code{n.analogs = 1}),
#'    otherwise it will give an error. Note that the analog dates are different for each member in case of 
#'    multimember downscaling, and are returned as a list, each element of the list corresponding to one member. 
#' 
#' @details 
#' 
#' \strong{Spatial consistency}
#' 
#' Several checks of spatial consistency are performed. In particular, note that both 'pred' (reanalysis) and 'sim' (model
#' simulations) should be in the same grid. This consistency must be ensured by the user prior to entering these arguments,
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
#' \strong{Construction of analogs using multiple neighbours}
#' 
#' The argument \code{sel.fun} controls how the analogs are constructed when considering more than the first neighbour (argument
#' \code{n.analogs} > 1). In this case the \code{"random"} choice randomly selects one of the \code{n.analogs} neighbours,
#'  while the \code{"mean"} choice will compute their average.
#' @seealso \code{\link{prinComp}} for details on principal component/EOF analysis
#' \code{\link{makeMultiGrid}} for multigrid creation
#' \code{loadGridData} and \code{loadStationData}, from package \pkg{loadeR}, for loading grids and station data respectively.
#' @importFrom fields rdist
#' @importFrom abind abind
#' @family downscaling
#' @references 
#' Benestad, R.E., Hanssen-Bauer, I. and Chen, D., 2008. Empirical-Statistical Downscaling,
#'  1st ed. World Scientific Publishing, Singapore
#'  
#' Gutierrez, J.M. \emph{et al.}, 2013. Reassessing Statistical downscaling techniques for
#'  their robust application under climate change conditions. J. Clim. 26, 171-188
#'  
#' Bedia, J. \emph{et al.}, 2013. Robust projections of Fire Weather Index in the Mediterranean
#'  using statistical downscaling. Clim. Change 120, 229-247.
#' 
#' @author J. Bedia 
#' @keywords internal

analogs <- function(y = y,
                    modelPars,
                    n.analogs = n.analogs,
                    sel.fun = sel.fun,
                    analog.dates = analog.dates,
                    parallel.pars = parallel.pars) {
      # modelPars <- ppModelSetup(y, pred, sim)
      #pred <- NULL
      if (isTRUE(modelPars$multi.member)) {
            y[[5]] <-  modelPars$init.dates
            y[[6]] <- modelPars$member.names
            names(y)[5:6] <- c("InitializationDates", "Members")
      } 
      #sim <- NULL
      n.analogs <- as.integer(n.analogs)
      if (n.analogs < 1) {
            stop("A minimum of 1 nearest analog must be selected in 'n.analogs'")
      }
      if (isTRUE(analog.dates) & n.analogs > 1) {
            stop("Analog dates are only returned for 1-neighbour analogs\n Set argument 'n.analogs = 1'")
      }
      sel.fun <- match.arg(sel.fun, choices = c("random", "mean"))
      # Analog search      
      message("[", Sys.time(), "] Calculating analogs ...")
      d.list <- if (parallel.pars$hasparallel) {
            on.exit(parallel::stopCluster(parallel.pars$cl))
            parallel::parLapply(cl = parallel.pars$cl, 1:length(modelPars$sim.mat), function(k) {
                  aux <- rdist(modelPars$sim.mat[[k]], modelPars$pred.mat)
                  apply(aux, 1, function(vec, n.analogs) {sort(vec, index.return = TRUE)$ix[1:n.analogs]}, n.analogs)
            })
      } else {
            lapply(1:length(modelPars$sim.mat), function(k) {
                  aux <- rdist(modelPars$sim.mat[[k]], modelPars$pred.mat)
                  apply(aux, 1, function(vec, n.analogs) {sort(vec, index.return = TRUE)$ix[1:n.analogs]}, n.analogs)
            })
      }
      modelPars$pred.mat <- NULL
      # Analog dates
      if (isTRUE(analog.dates)) {
            analog.date.list <- lapply(1:length(d.list), function(k) y$Dates$start[d.list[[k]]])
            attr(y, "analog.dates") <- analog.date.list
            analog.date.list <- NULL
      }
      # Analog assignation
      if (isTRUE(modelPars$stations)) {
            if (n.analogs > 1) {
                  out.list <- lapply(1:length(d.list), function(k) {
                        aux.mat <- matrix(nrow = ncol(d.list[[k]]), ncol = dim(y$Data)[grep("station", attr(y$Data, "dimensions"))])
                        for (i in 1:nrow(aux.mat)) {
                              aux <- y$Data[d.list[[k]][ ,i], ]
                              aux.mat[i, ] <- switch(sel.fun, 
                                                     "random" = aux[sample(1:nrow(aux), 1), ],
                                                     "mean" = apply(aux, 2, mean, na.rm = TRUE))
                        }
                        return(aux.mat)
                  })
                  message("[", Sys.time(), "] Done.")   
            } else {
                  out.list <- lapply(1:length(d.list), function(k) y$Data[d.list[[k]], ])
                  message("[", Sys.time(), "] Done.")   
            }
      } else {
            obs.mat <- array3Dto2Dmat(y$Data)
            obs.coords <- getCoordinates(y)
            if (n.analogs > 1) {
                  out.list <- lapply(1:length(d.list), function(k) {
                        aux.mat <- matrix(nrow = ncol(d.list[[k]]), ncol = ncol(obs.mat))
                        for (i in 1:ncol(d.list[[k]])) {
                              aux <- obs.mat[d.list[[k]][ ,i], ]
                              aux.mat[i, ] <- switch(sel.fun, 
                                                     "random" = aux[sample(1:nrow(aux), 1), ],
                                                     "mean" = apply(aux, 2, mean, na.rm = TRUE))
                        }
                        aux.mat <- mat2Dto3Darray(aux.mat, obs.coords$x, obs.coords$y)
                        return(aux.mat)
                  })
                  message("[", Sys.time(), "] Done.")   
            } else {
                  out.list <- lapply(1:length(d.list), function(k) {
                        mat2Dto3Darray(obs.mat[d.list[[k]], ], obs.coords$x, obs.coords$y)
                  })
                  message("[", Sys.time(), "] Done.")   
            }
            obs.mat <- NULL
      }
      d.list <- NULL
      # Data array
      #obs$Data 
      o <- drop(unname(do.call("abind", c(out.list, along = -1))))
      out.list <- NULL
      # New data attributes
      # Date replacement
      return(o)
}
# End

