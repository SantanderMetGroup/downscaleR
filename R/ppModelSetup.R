#' @title Data preparation for perfect-prog model construction
#' @description Performs the various data pre-processing and formatting steps required
#' to implement the different perfect-prog downscaling methods.
#' @template templateObsPredSim
#' @return A list with differet components required by the downscaling functions:
#' \itemize{
#' \item \code{stations} Logical indicating whether the predictand comes from station data
#' \item \code{multi.member} Logical indicating whether the simulation data is multi-member or not
# #' \item \code{obs} Same as input parameter of the same name
#' \item \code{pred.mat} 2D matrix of predictors. Predictors are arranged in columns and time in rows
#' \item \code{sim.mat} Either a single 2D matrix or a list of them in case of multi-member predictions,
#' containing the simulation data to compute the downscaling predictions.
#' \item \code{sim.dates} The Dates element of the simulation data (either a list of start/end lists for
#' several predictors or a list of two with the estart/end of a single predictor. See \code{\link{dateReplacement}}
#' for details).
#' }
#' @details The function accepts either PCA or raw fields as predictors. In the first case, it handles the 
#' mapping of EOFs onto the simulation fields.
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal
#' @export
#' @importFrom abind asub

ppModelSetUp <- function(obs, pred, sim) {
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
                  multi.member <- TRUE
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
                  multi.member <- FALSE
                  simsc.list <- lapply(1:n.vars, function(x) {
                        attr(simsc.list[[x]], "dimensions") <- dimNames.sim[-var.dim.index]
                        aux <- array3Dto2Dmat(simsc.list[[x]])
                        aux <- (aux - mu.list[[x]]) / sigma.list[[x]]
                        return(aux)
                  })
            }
      } else { # Field
            if ("member" %in% dimNames.sim) { # Multimember field
                  multi.member <- TRUE
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
                  multi.member <- FALSE
                  aux <- array3Dto2Dmat(sim$Data)
                  aux <- (aux - mu.list[[1]]) / sigma.list[[1]]
                  simsc.list <- list(aux)
                  aux <- NULL
            }
      }
      # Sim 2D matrix
      if (length(simsc.list) > 1) {
            if ("member" %in% attr(sim$Data, "dimensions")) {
                  multi.member <- TRUE
                  sim.mat <- rep(list(bquote()), n.mem)
                  for (i in 1:n.mem) {
                        aux <- lapply(1:length(simsc.list), function(x) simsc.list[[x]][[i]])
                        sim.mat[[i]] <- do.call("cbind", aux)
                  }
            } else {
                  multi.member <- FALSE
                  sim.mat <- list(do.call("cbind", simsc.list))
            }
      } else {
            sim.mat <- simsc.list
      }
      simsc.list <- NULL
      attr(sim.mat, "dataset") <- attributes(sim)$"dataset"
      sim.dates <- sim$Dates
      sim <- NULL
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
      return(list("stations" = stations, "multi.member" = multi.member, "pred.mat" = pred.mat, "sim.mat" = sim.mat, "sim.dates" = sim.dates))
}
# End      
