#' @title Analog downscaling
#' 
#' @description Implementation of the downscaling analogs method
#' 
#' @template templateObsPredSim
#' @param n.neigh Integer indicating the number of closest neigbours to retain for analog construction. Default to 1.
#' @param sel.fun Criterion for the construction of analogs when several neigbours are chosen. Ignored when \code{n.neig = 1}.
#' Current values are \code{"random"} (the default) and \code{"mean"}. See details.
#' 
#' @details 
#' 
#' \strong{Spatial consistency}
#' 
#' Several checks of spatial consistency are performed. In particular, note that both 'pred' (reanalysis) and 'sim' (model
#' simulations) should be in the same grid. This consistency must be ensured by the user prior to entering these arguments,
#' for instance by means of the \code{\link{interpGridData}} function in conjunction with the \code{\link{getGrid}} method.
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
#' \code{n.neigh} > 1). In this case the \code{"random"} choice randomly selects one of the \code{n.neigh} neighbours,
#'  while the \code{"mean"} choice will compute their average.
#' 
#' @seealso \code{\link{prinComp}} for details on principal component/EOF analysis
#' \code{\link{loadMultiField}}, \code{\link{makeMultiField}} for multifield creation
#' \code{\link{loadGridData}} and \code{\link{loadStationData}} for loading fields and station data respectively.
#' 
#' @export
#' 
#' @importFrom fields rdist
#' @importFrom abind abind
#' 
#' @family downscaling
#' 
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
#' @author J. Bedia \email{joaquin.bedia@@gmail.com} 
#'

analogs <- function(obs, pred, sim, n.neigh = 1, sel.fun = c("random", "mean")) {
      modelPars <- ppModelSetUp(obs, pred, sim)
      pred <- NULL
      sim <- NULL
      n.neigh <- as.integer(n.neigh)
      if (n.neigh < 1) {
            stop("A minimum of 1 nearest neighbour must be selected in 'n.neigh'")
      }
      sel.fun <- match.arg(sel.fun, choices = c("random", "mean"))
      # Analog search      
      message("[", Sys.time(), "] Calculating analogs ...")
      d.list <- lapply(1:length(modelPars$sim.mat), function(x) {
            aux <- rdist(modelPars$sim.mat[[x]], modelPars$pred.mat)
            aux <- apply(aux, 1, function(vec, n.neigh) {sort(vec, index.return = TRUE)$ix[1:n.neigh]}, n.neigh)
            return(aux)
      })
      modelPars$pred.mat <- NULL
      # Analog assignation
      if (isTRUE(modelPars$stations)) {
            if (n.neigh > 1) {
                  out.list <- lapply(1:length(d.list), function(x) {
                        aux.mat <- matrix(nrow = ncol(d.list[[x]]), ncol = dim(obs$Data)[grep("station", attr(obs$Data, "dimensions"))])
                        for (i in 1:nrow(aux.mat)) {
                              aux <- obs$Data[d.list[[x]][ ,i], ]
                              aux.mat[i, ] <- switch(sel.fun, 
                                                     "random" = aux[sample(1:nrow(aux), 1), ],
                                                     "mean" = apply(aux, 2, mean, na.rm = TRUE))
                        }
                        return(aux.mat)
                  })
                  message("[", Sys.time(), "] Done.")   
            } else {
                  out.list <- lapply(1:length(d.list), function(x) {obs$Data[d.list[[x]], ]})
                  message("[", Sys.time(), "] Done.")   
            }
      } else {
            obs.mat <- array3Dto2Dmat(obs$Data)
            if (n.neigh > 1) {
                  out.list <- lapply(1:length(d.list), function(x) {
                        aux.mat <- matrix(nrow = ncol(d.list[[x]]), ncol = ncol(obs.mat))
                        for (i in 1:ncol(d.list[[x]])) {
                              aux <- obs.mat[d.list[[x]][ ,i], ]
                              aux.mat[i, ] <- switch(sel.fun, 
                                                     "random" = aux[sample(1:nrow(aux), 1), ],
                                                     "mean" = apply(aux, 2, mean, na.rm = TRUE))
                        }
                        aux.mat <- mat2Dto3Darray(aux.mat, x.obs, y.obs)
                        return(aux.mat)
                  })
                  message("[", Sys.time(), "] Done.")   
            } else {
                  out.list <- lapply(1:length(d.list), function(x) {
                        mat2Dto3Darray(obs.mat[d.list[[x]], ], x.obs, y.obs)
                  })
                  message("[", Sys.time(), "] Done.")   
            }
            obs.mat <- NULL
      }
      d.list <- NULL
      # Data array
      dimNames <- attr(obs$Data, "dimensions")
      # Remove "station" from dimensions for single-station objects
      st.dim.index <- grep("station", dimNames)
      if (!identical(st.dim.index, integer(0))) {
            dim.st <- dim(obs$Data)[st.dim.index]
            if (identical(dim.st, 1L)) {
                  dimNames <- dimNames[-st.dim.index]
            }
      }
      obs$Data <- drop(unname(do.call("abind", c(out.list, along = -1))))
      out.list <- NULL
      # if ("member" %in% attr(sim$Data, "dimensions")) {
      if (isTRUE(modelPars$multi.member)) {
            attr(obs$Data, "dimensions") <- c("member", dimNames)
      } else {
            attr(obs$Data, "dimensions") <- dimNames
      }
      # New data attributes      
      attr(obs$Data, "downscaling:method") <- "analogs"
      attr(obs$Data, "downscaling:simulation_data") <- attr(modelPars$sim.mat, "dataset")
      # Date replacement
      obs$Dates <- dateReplacement(obs$Dates, modelPars$sim.dates) 
      return(obs)
}
# End
# End

