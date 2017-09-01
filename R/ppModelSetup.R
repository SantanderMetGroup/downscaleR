#' @title Data preparation for perfect-prog model construction
#' @description Performs the various data pre-processing and formatting steps required
#' to implement the different perfect-prog downscaling methods.
#' @template templateObsPredSim
#' @return A list with differet components required by the downscaling functions:
#' \itemize{
#' \item \code{stations} Logical indicating whether the predictand comes from station data
#' \item \code{multi.member} Logical indicating whether newdata is multi-member or not
#  \item \code{obs} Same as input parameter of the same name
#' \item \code{pred.mat} 2D matrix of predictors. Predictors are arranged in columns and time in rows
#' \item \code{sim.mat} A list of 2D matrices containing the prediction data. The list is of length \emph{n},
#' being \emph{n} the number of members considered (n = 1 for deterministic/single member predictions).
#' \item \code{sim.dates} The Dates element of the newdata (either a list of start/end lists for
#' several predictors or a list of two with the estart/end of a single predictor. See \code{\link{dateReplacement}}
#' for details).
#' \item \code{init.dates} Initialization dates inherited from 'newdata' if a forecast. NULL otherwise.
#' \item \code{member.names} Names of the members inherited from 'newdata' if a forecast. NULL otherwise.
#' }
#' @details The function accepts either PCA or raw grids as predictors. In the first case, it handles the 
#' mapping of EOFs onto the newdata grids.
#' @author J. Bedia 
#' @keywords internal
#' @importFrom abind asub adrop
#' @importFrom stats sd
#' @importFrom transformeR array3Dto2Dmat getDim getShape
#' @importFrom magrittr %>%

ppModelSetup <- function(y, x, newdata) {
      stations <- ifelse("loc" %in% getDim(y), TRUE, FALSE)
      if ("scaled:method" %in% names(attributes(x))) {
            use.PCs <- TRUE
            x.pred <- attr(x, "xCoords")
            y.pred <- attr(x, "yCoords")
            time.pred <- attr(x, "dates_start")
      } else {
            if (length(getShape(x)) < 3) {
                  stop("'x' must be a grid/multigrid\nSingle point selections are not allowed", call. = FALSE)
            }
            use.PCs <- FALSE
            x.pred <- getCoordinates(x)$x
            y.pred <- getCoordinates(x)$y
            time.pred <- getRefDates(x, "start")
      }
      # Georef simulations
      if (length(getShape(newdata)) < 3) {
          stop("'newdata' must be a grid/multigrid\nSingle point selections are not allowed", call. = FALSE)
      }
      x.sim <- newdata$xyCoords$x
      y.sim <- newdata$xyCoords$y
      # Spatial consistency check x-newdata
      if (!isTRUE(all.equal(x.pred, x.sim, tolerance = 1e-03)) || !isTRUE(all.equal(y.pred, y.sim, tolerance = 1e-03))) {
          stop("'x' and 'newdata' datasets are not spatially consistent", call. = FALSE)
      }
      # Temporal matching check (y-x)
      if (!identical(as.POSIXlt(y$Dates$start)$yday, as.POSIXlt(time.pred)$yday) || !identical(as.POSIXlt(y$Dates$start)$year, as.POSIXlt(time.pred)$year)) {
            stop("Observed and predicted time series should match in start/end and length", call. = FALSE)
      }     
      # Date replacement
      new.dates <- dateReplacement(y$Dates, newdata$Dates)
      y <- time.pred <- NULL
      # Number of variables x-newdata
      n.vars <- if (isTRUE(use.PCs)) {
            ifelse(length(x) > 1, length(x) - 1, 1)
      } else {
            length(x$Variable$varName)
      }
      if (length(newdata$Variable$varName) != n.vars) {
            stop("The number of variables of predictor (x) and newdata does not match", call. = FALSE)      
      }
      # Scaling and centering of simulation data. Parameters are inherited from the predictors
      if (isTRUE(use.PCs)) {
            mu.list <- lapply(1:n.vars, function(k) attributes(x[[k]][[1]]$orig)$"scaled:center")
            sigma.list <- lapply(1:n.vars, function(k) attributes(x[[k]][[1]]$orig)$"scaled:scale")
            pred.mat <- x$COMBINED[[1]]$PCs
            if (is.null(pred.mat)) pred.mat <- x[[1]][[1]]$PCs
      } else {
            x <- redim(x, member = FALSE, var = TRUE)
            # x rescaled matrix
            n.vars <- getShape(x, "var")
            varNames <- x$Variable$varName
            Xsc.list <- lapply(1:n.vars, function(idx) {
                  aux <- suppressWarnings(subsetGrid(x, var = varNames[idx], drop = TRUE)[["Data"]]) %>% array3Dto2Dmat()
                  mu <- mean(aux, na.rm = TRUE)
                  sigma <- sd(aux, na.rm = TRUE)
                  aux <- (aux - mu) / sigma
                  attr(aux, "scaled:center") <- mu
                  attr(aux, "scaled:scale") <- sigma
                  return(aux)
            })
            # Scaling parameters
            mu.list <- lapply(1:n.vars, function(k) attributes(Xsc.list[[k]])$"scaled:center")
            sigma.list <- lapply(1:n.vars, function(k) attributes(Xsc.list[[k]])$"scaled:scale")
            # standardized x 2D matrix
            pred.mat <- do.call("cbind", Xsc.list)
            Xsc.list <- NULL
      }
      # Scaling and centering of simulation data
      newdata <- redim(newdata, var = TRUE, member = TRUE)
      n.mem <- getShape(newdata, "member")
      varNames <- newdata$Variable$varName
      simsc.list <- lapply(1:n.vars, function(idx) {
            a <- suppressWarnings(subsetGrid(newdata, var = varNames[idx], drop = TRUE)) %>% redim(member = TRUE)
            lapply(1:n.mem, function(id.mem) {
                  aux <- subsetGrid(a, members = id.mem, drop = TRUE)[["Data"]] %>% array3Dto2Dmat()
                  (aux - mu.list[[idx]]) / sigma.list[[idx]]
            })
      })
      # newdata 2D matrix
      sim.mat <- vector("list", n.mem)
      for (i in 1:n.mem) { 
            aux <- lapply(1:length(simsc.list), function(k) simsc.list[[k]][[i]])
            sim.mat[[i]] <- do.call("cbind", aux)
      }
      simsc.list <- NULL
      sim.dataset <- attr(newdata, "dataset")
      init.dates <- mems <- NULL
      if (!is.null(newdata$InitializationDates)) init.dates <- newdata$InitializationDates
      if (!is.null(newdata$Members)) mems <- newdata$Members
      newdata <- NULL
      # Projection of simulated grid onto the predictor EOFs       
      if (isTRUE(use.PCs)) {
            sim.mat <- tryCatch(expr = lapply(1:length(sim.mat), function(k) {
                  sim.mat[[k]] %*% x$COMBINED[[1]]$EOFs
            }), error = function(err) {
                  lapply(1:length(sim.mat), function(k) {
                        sim.mat[[k]] %*% x[[1]][[1]]$EOFs 
                  })
            })
      }
      x <- NULL
      multi.member <- ifelse(is.null(mems), FALSE, TRUE)
      return(list("stations" = stations,
                  "multi.member" = multi.member,
                  "pred.mat" = pred.mat,
                  "sim.mat" = sim.mat,
                  "sim.dates" = new.dates,
                  "sim.dataset" = sim.dataset,
                  "init.dates" = init.dates,
                  "member.names" = mems)
      )
}
# End      
