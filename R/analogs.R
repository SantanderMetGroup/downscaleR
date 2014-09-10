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
#' @importFrom abind asub
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
      n.neigh <- as.integer(n.neigh)
      if (n.neigh < 1) {
            stop("A minimum of 1 nearest neighbour must be selected in 'n.neigh'")
      }
      sel.fun <- match.arg(sel.fun, choices = c("random", "mean"))
      # Georef predictand
      if ("station" %in% attr(obs$Data, "dimensions")) {
            stations <- TRUE
            x.obs <- obs$xyCoords[ ,1]
            y.obs <- obs$xyCoords[ ,2]
      } else {
            stations <- FALSE
            x.obs <- obs$xyCoords$x
            y.obs <- obs$xyCoords$y
      }
      # Georef predictors
      if ("scaled:method" %in% names(attributes(pred))) {
            use.PCs <- TRUE
            x.pred <- attr(pred, "xCoords")
            y.pred <- attr(pred, "yCoords")
            time.pred <- attr(pred, "dates_start")
      } else {
            if (length(dim(pred$Data)) < 3) {
                  stop("'pred' must be a field/multifield\nSingle point selections are not allowed")
            }
            use.PCs <- FALSE
            x.pred <- pred$xyCoords$x
            y.pred <- pred$xyCoords$y
            if (is.null(names(pred$Dates))) {
                  time.pred <- pred$Dates[[1]]$start   
            } else {
                  time.pred <- pred$Dates$start 
            }
      }
      # Georef simulations
      if (length(dim(sim$Data)) < 3) {
            stop("'pred' must be a field/multifield\nSingle point selections are not allowed")
      }
      x.sim <- sim$xyCoords$x
      y.sim <- sim$xyCoords$y
      # Spatial consistency check pred-sim
      if (!isTRUE(all.equal(x.pred, x.sim, tolerance = 1e-03)) || !isTRUE(all.equal(y.pred, y.sim, tolerance = 1e-03))) {
            stop("'pred' and 'sim' datasets are not spatially consistent")
      }
      # Temporal matching check (obs-pred)
      if (!identical(as.POSIXlt(obs$Dates$start)$yday, as.POSIXlt(time.pred)$yday) || !identical(as.POSIXlt(obs$Dates$start)$year, as.POSIXlt(time.pred)$year)) {
            stop("Observed and predicted time series should match in start/end and length")
      }     
      time.pred <- NULL
      # Number of variables pred-sim
      if (isTRUE(use.PCs)) {
            n.vars <- ifelse(length(pred) > 1, length(pred) - 1, 1)
      } else {
            n.vars <- length(pred$Variable$varName)
      }
      if (length(sim$Variable$varName) != n.vars) {
            stop("The number of variables of predictor and simulation datasets does not match")      
      }
      # Index of the closest pred/sim points to the observations
#       lon.index <- unlist(lapply(1:length(x.obs), function(x) which.min((x.obs[x] - x.pred)^2)))
#       lat.index <- unlist(lapply(1:length(y.obs), function(x) which.min((y.obs[x] - y.pred)^2)))
      # Scaling and centering of simulation data. Parameters are inherited from the predictors
      if (isTRUE(use.PCs)) {
            mu.list <- lapply(1:n.vars, function(x) {attributes(pred[[x]][[1]])$"scaled:center"})
            sigma.list <- lapply(1:n.vars, function(x) {attributes(pred[[x]][[1]])$"scaled:scale"})
            pred.mat <- pred$COMBINED[[1]]$PCs
            if (is.null(pred.mat)) {
                  pred.mat <- pred[[1]][[1]]$PCs
            }      
      } else {
            # pred rescaled matrix
            dimNames <- attr(pred$Data, "dimensions")
            if ("var" %in% attr(pred$Data, "dimensions")) { # multifield
                  var.dim.index <- grep("var", dimNames)
                  n.vars <- dim(pred$Data)[var.dim.index]
                  Xsc.list <- lapply(1:n.vars, function(idx) {
                        aux <- asub(pred$Data, idx, var.dim.index)
                        attr(aux, "dimensions") <- dimNames[-var.dim.index]
                        aux <- array3Dto2Dmat(aux)   
                        mu <- mean(aux, na.rm = TRUE)
                        sigma <- sd(aux, na.rm = TRUE)
                        aux <- (aux - mu) / sigma
                        attr(aux, "scaled:center") <- mu
                        attr(aux, "scaled:scale") <- sigma
                        return(aux)
                  })
            } else {
                  n.vars <- 1
                  aux <- array3Dto2Dmat(pred$Data)
                  mu <- mean(aux, na.rm = TRUE)
                  sigma <- sd(aux, na.rm = TRUE)
                  aux <- (aux - mu) / sigma
                  attr(aux, "scaled:center") <- mu
                  attr(aux, "scaled:scale") <- sigma
                  Xsc.list <- list(aux)
                  aux <- NULL
            }
            # Scaling parameters
            mu.list <- lapply(1:n.vars, function(x) attributes(Xsc.list[[x]])$"scaled:center")
            sigma.list <- lapply(1:n.vars, function(x) attributes(Xsc.list[[x]])$"scaled:scale")
            # standardized pred 2D matrix
            pred.mat <- do.call("cbind", Xsc.list)
            Xsc.list <- NULL
      }
      # Scaling and centering of simulation data
      dimNames.sim <- attr(sim$Data, "dimensions")
      if ("var" %in% dimNames.sim) { # Multifield
            var.dim.index <- grep("var", dimNames.sim)
            simsc.list <- lapply(1:n.vars, function(idx) {asub(sim$Data, idx, var.dim.index)})
            if ("member" %in% dimNames.sim) { # Multifield multimember
                  mem.dim.index <- grep("member", dimNames.sim[-var.dim.index])
                  n.mem <- dim(sim$Data)[-var.dim.index][mem.dim.index]
                  for (i in 1:n.vars) {
                        simsc.list[[i]] <- lapply(1:n.mem, function(id.mem) {
                              aux <- asub(simsc.list[[i]], id.mem, mem.dim.index)
                              attr(aux, "dimensions") <- dimNames.sim[-match(c("var","member"), dimNames.sim)]
                              aux <- array3Dto2Dmat(aux)
                              aux <- (aux - mu.list[[i]]) / sigma.list[[i]]
                              return(aux)
                        })
                  }
            } else { # Multifield (no members)
                  simsc.list <- lapply(1:n.vars, function(x) {
                        attr(simsc.list[[x]], "dimensions") <- dimNames.sim[-var.dim.index]
                        aux <- array3Dto2Dmat(simsc.list[[x]])
                        aux <- (aux - mu.list[[x]]) / sigma.list[[x]]
                        return(aux)
                  })
            }
      } else { # Field
            if ("member" %in% dimNames.sim) { # Multimember field
                  mem.dim.index <- grep("member", dimNames.sim)
                  n.mem <- dim(sim$Data)[mem.dim.index]
                  simsc.list <- lapply(1:length(n.mem), function(x) {
                        aux <- asub(sim$Data, n.mem, mem.dim.index)
                        attr(aux, "dimensions") <- dimNames.sim[-mem.dim.index]
                        aux <- array3Dto2Dmat(aux)
                        aux <- (aux - mu.list[[x]]) / sigma.list[[x]]
                        return(aux)
                  })
            } else { # Field (no multimember)
                  aux <- array3Dto2Dmat(sim$Data)
                  aux <- (aux - mu.list[[1]]) / sigma.list[[1]]
                  simsc.list <- list(aux)
                  aux <- NULL
            }
      }
      # Sim 2D matrix
      if (length(simsc.list) > 1) {
            if ("member" %in% attr(sim$Data, "dimensions")) {
                  sim.mat <- rep(list(bquote()), n.mem)
                  for (i in 1:n.mem) {
                        aux <- lapply(1:length(simsc.list), function(x) simsc.list[[x]][[i]])
                        sim.mat[[i]] <- do.call("cbind", aux)
                  }
            } else {
                  sim.mat <- list(do.call("cbind", simsc.list))
            }
      } else {
            sim.mat <- simsc.list
      }
      simsc.list <- NULL
      # Projection of simulated field onto the predictor EOFs       
      if (isTRUE(use.PCs)) {
            sim.mat <- tryCatch(expr = lapply(1:length(sim.mat), function(x) {
                        t(t(pred$COMBINED[[1]]$EOFs) %*% t(sim.mat[[x]]))
                  }), error = function(err) {
                        l <- lapply(1:length(sim.mat), function(x) {t(t(pred[[1]][[1]]$EOFs) %*% t(sim.mat[[x]]))})
                        return(l)
                  })
      }
      pred <- NULL
      # Analog search      
      message("[", Sys.time(), "] Calculating analogs ...")
      d.list <- lapply(1:length(sim.mat), function(x) {
            aux <- rdist(sim.mat[[x]], pred.mat)
            aux <- apply(aux, 1, function(vec, n.neigh) {sort(vec, index.return = TRUE)$ix[1:n.neigh]}, n.neigh)
            return(aux)
      })
      pred.mat <- NULL
      # Analog assignation
      if (isTRUE(stations)) {
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
      if ("member" %in% attr(sim$Data, "dimensions")) {
            attr(obs$Data, "dimensions") <- c("member", dimNames)
      } else {
            attr(obs$Data, "dimensions") <- dimNames
      }
      # New data attributes      
      attr(obs$Data, "downscaling:method") <- "analogs"
      attr(obs$Data, "downscaling:simulation_data") <- attributes(sim)$"dataset"
      # Date replacement
      obs$Dates <- dateReplacement(obs$Dates, sim$Dates) 
      return(obs)
}
# End

