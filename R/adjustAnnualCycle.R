#' @title Annual cycle scaling of simulation data
#' 
#' @description Annual cycle scaling of simulation data w.r.t. the predictors
#' 
#' @param pred Predictor field
#' @param sim Simulation field
#' @param ref Optional reference field. A field from where the scaling and centering parameters are taken. See details.
#' @param ensemble Logical flag. Should the correction of the mean be performed w.r.t. the ensemble mean (\code{TRUE})
#'  or member by member independently (\code{ensemble = FALSE})?. Ignored if \code{ref} (and thus \code{sim}) are not
#'   multimembers. Default to \code{FALSE}.
#'  
#' @return A field with the rescaled simulation data (\code{sim}), with the centering parameters
#'  indicated as an attribute
#'  
#' @details The reference field is used to correct the simulation (test) data, as follows:
#' 
#' \deqn{sim' = sim - mu_ref + mu_pred}
#' 
#' , where \emph{mu} corresponds to the monthly climatological mean considering the training period,
#' and \emph{sim'} is the corrected simulation (test) data. The way \emph{mu_ref} is computed in case
#' of multimember fields is controlled by the argument \code{ensemble}.
#' 
#' The \code{ref} usually corresponds to the control run of the GCM in the training period in climate change applications,
#' or the hindcast data for the training period in s2d applications. Note that by default \code{ref = NULL}. In this 
#' case it will be assumed to be the \code{pred} field. This can be used for instance when train and test correspond
#' to the same model.
#' 
#' @keywords internal
#' @export
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}



adjustAnnualCycle <- function(pred, sim, ref = NULL, ensemble = FALSE) {
      use.ref <- ifelse(is.null(ref), FALSE, TRUE)
      if (is.null(ref)) ref <- pred
      dimNames <- attr(ref$Data, "dimensions")
      if(!identical(dimNames, attr(sim$Data, "dimensions"))) stop("Input and reference field dimensions do not match")
      seas <- getSeason(pred)
      if (!identical(seas, getSeason(ref))) stop("Season of input and reference fields do not match")
      var.names <- ref$Variable$varName
      if (!identical(var.names, sim$Variable$varName) | !identical(var.names, pred$Variable$varName)) stop("Variable(s) of predictor and simulation fields do not match")
      aux.ind <- grep(paste(c("var","member","lat","lon"), collapse = "|"), dimNames)
      if (!identical(dim(ref$Data)[aux.ind], dim(sim$Data)[aux.ind])) stop("Spatial and/or ensemble dimensions of sim and reference fields do not match")
      mon <- if ("var" %in% dimNames) {
            as.POSIXlt(sim$Dates[[1]]$start)$mon + 1
      } else {
            as.POSIXlt(sim$Dates$start)$mon + 1
      }
      index <- unlist(sapply(seas, function(x) which(mon == x), simplify = TRUE))
      mem.ind <- grep("member", attr(ref$Data, "dimensions"))
      n.mem <- ifelse(length(mem.ind) > 0, dim(ref$Data)[mem.ind], 1L)
      message("[", Sys.time(), "] Calculating parameters...")
      # MU: monthly mean pars of the predictor
      center.list.pred <- lapply(1:length(var.names), function(x) {
            a <- suppressWarnings(subsetField(pred, var = var.names[x]))
            b <- lapply(seas, function(y) {
                  aux <- subsetField(a, season = y)$Data
                  aux <- apply(aux, grep("member", attr(aux, "dimensions"), invert = TRUE), mean, na.rm = TRUE)
                  attr(aux, "dimensions") <- c("time", "lat", "lon")                                    
                  colMeans(array3Dto2Dmat(aux))
            })
            names(b) <- month.abb[seas]
            return(b)
      })
      names(center.list.pred) <- var.names
      center.list.pred <- rep(list(center.list.pred), n.mem)
      # MU': monthly mean pars of the reference field
      if (!use.ref) {
            center.list.ref <- center.list.pred
      } else {
            if (n.mem > 1) {
                  if (isTRUE(ensemble)) {
                        aux1 <- lapply(1:length(var.names), function(x) {
                              a <- suppressWarnings(subsetField(ref, var = var.names[x]))
                              b <- lapply(seas, function(y) {
                                    aux <- subsetField(a, season = y)$Data
                                    aux <- apply(aux, grep("member", attr(aux, "dimensions"), invert = TRUE), mean, na.rm = TRUE)
                                    attr(aux, "dimensions") <- c("time", "lat", "lon")                                    
                                    colMeans(array3Dto2Dmat(aux))
                              })
                              names(b) <- month.abb[seas]
                              return(b)
                        })
                        names(aux1) <- var.names
                        center.list.ref <- rep(list(aux1), n.mem)
                  } else {
                        center.list.ref <- lapply(1:n.mem, function(x) {
                              a <- subsetField(ref, members = x)
                              aux1 <- lapply(1:length(var.names), function(y) {
                                    b <- suppressWarnings(subsetField(a, var = var.names[y]))
                                    aux <- lapply(seas, function(z) {
                                          colMeans(array3Dto2Dmat(subsetField(b, season = z)$Data))
                                    })
                                    names(aux) <- month.abb[seas]
                                    return(aux)
                              })
                              names(aux1) <- var.names
                              return(aux1)
                        })
                  }
                  names(center.list.ref) <- ref$Members
            } else {
                  center.list.ref <- lapply(1:length(var.names), function(x) {
                        a <- suppressWarnings(subsetField(ref, var = var.names[x]))
                        b <- lapply(seas, function(y) {
                              aux <- subsetField(a, season = y)$Data
                              aux <- apply(aux, grep("member", attr(aux, "dimensions"), invert = TRUE), mean, na.rm = TRUE)
                              attr(aux, "dimensions") <- c("time", "lat", "lon")                                    
                              colMeans(array3Dto2Dmat(aux))
                        })
                        names(b) <- month.abb[seas]
                        return(b)
                  })
                  names(center.list.ref) <- var.names
                  center.list.ref <- list(center.list.ref)
            }
      }
      ref <- pred <- NULL
      # In base::scale, each column of x has the corresponding value from center **subtracted** from it.
      # Thus: center = MUpred - MUref
      message("[", Sys.time(), "] Rescaling...")
      aux.list <- lapply(1:length(var.names), function(v) {
            sf.var <- suppressWarnings(subsetField(sim, var = var.names[v]))
            l <- lapply(1:n.mem, function(m) {
                  aux <- array3Dto2Dmat(suppressWarnings(subsetField(sf.var, members = m))$Data)
                  mu.list <- lapply(1:length(seas), function(s) {
                        scale(aux[which(mon == seas[s]), ], 
#                               center = (center.list.ref[[m]][[v]][[s]] + center.list.pred[[m]][[v]][[s]]),
                              center = (center.list.ref[[m]][[v]][[s]] - center.list.pred[[m]][[v]][[s]]),
                              scale = FALSE)
                  })
                  arr3d <- mat2Dto3Darray(do.call("rbind", mu.list)[index, ], sim$xyCoords$x, sim$xyCoords$y)
                  ll <- sapply(1:length(mu.list), function(i) attributes(mu.list[[i]])[-1])
                  names(ll) <- month.name[seas]
                  return(list(ll, arr3d))
            })
            arr4d <- unname(do.call("abind", c(lapply(1:length(l), function(x) l[[x]][[2]]), along = -1L)))
            par.list <- lapply(1:length(l), function(x) l[[x]][[1]])
            return(list(par.list, arr4d))
      })
      message("[", Sys.time(), "] Done.")
      arr <- if (length(aux.list) > 1) {# Multifield
                  unname(do.call("abind", c(lapply(1:length(aux.list), function(x) aux.list[[x]][[2]]), along = -1L)))
      } else {
            aux.list[[1]][[2]]
      }
      par.list <- lapply(1:length(aux.list), function(x) aux.list[[x]][[1]])
      names(par.list) <- var.names
      attr(arr, "dimensions") <- dimNames
      attr(arr, "scale::center") <- par.list
      sim$Data <- arr
      return(sim)
}
# End 