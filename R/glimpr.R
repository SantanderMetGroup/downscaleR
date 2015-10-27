#' @title GLM downscaling
#' @description Implementation of GLM's to downscale precipitation data
#' @template templateObsPredSim
#' @param pr.threshold Value below which precipitation amount is considered zero 
#' @param n.eofs Integer indicating the number of EOFs to be used as predictors 
#' @family downscaling
#' @author J Bedia \email{joaquin.bedia@@gmail.com}
#' @importFrom abind abind

glimpr<- function(obs, pred, sim, pr.threshold = .1, n.eofs = 15) {
      modelPars <- ppModelSetup(obs, pred, sim)
      if(is.null(n.eofs)){
            n.eofs <-  dim(modelPars$pred.mat)[2]
      }
      pred <- NULL
      sim <- NULL
      # Prepare predictand matrix
      ymat <- if (isTRUE(modelPars$stations)) {
            obs$Data
      } else {
            if (is.null(dim(obs$Data))) {
                  as.matrix(obs$Data)
            } else {
                  array3Dto2Dmat(obs$Data)
            }
      }
      # binarize      
      ymat.bin <- ymat
      ymat.bin[which(ymat < pr.threshold)] <- 0
      ymat.bin[which(ymat > 0)] <- 1L
      wet.prob <- apply(ymat.bin, 2, function(x) {sum(x, na.rm = TRUE) / length(na.exclude(x))})
      # Filter empty gridcell series
      rm.ind <- which(is.na(wet.prob))
      # Valid columns (cases)
      cases <- if (length(rm.ind) > 0) {
            (1:ncol(ymat))[-rm.ind]      
      } else {
            1:ncol(ymat)
      }
      # Models
      message("[", Sys.time(), "] Fitting models ...")
      pred.list <- suppressWarnings(lapply(1:length(cases), function(x) {
            mod.bin <- glm(ymat.bin[ ,cases[x]] ~ ., data = as.data.frame(modelPars$pred.mat)[,1:n.eofs], family = binomial(link = "logit"))
            sims.bin <- sapply(1:length(modelPars$sim.mat), function(i) {
                
                  pred <- predict(mod.bin, newdata = as.data.frame(modelPars$sim.mat[[i]])[,1:n.eofs], type = "response")
                  pred[which(pred >= quantile(pred, 1 - wet.prob[cases[x]]))] <- 1L
                  pred <- as.integer(pred)
            })
            mod.bin <- NULL
            wet <- which(ymat[ ,cases[x]] != 0)
            # TODO: handle exception when model is over-specified (https://stat.ethz.ch/pipermail/r-help/2000-January/009833.html)
            mod.gamma <- glm(ymat[ ,cases[x]] ~., data = as.data.frame(modelPars$pred.mat)[,1:n.eofs], family = Gamma(link = "log"), subset = wet)
            sims <- sapply(1:length(modelPars$sim.mat), function(i) {
                  unname(predict(mod.gamma, newdata = as.data.frame(modelPars$sim.mat[[i]])[,1:n.eofs], type = "response") * sims.bin[ ,i])
            })
            return(unname(sims))
      }))
      # Refill with NaN series
      if (length(rm.ind) > 0) {
            aux.list <- rep(list(matrix(nrow = nrow(ymat), ncol = length(modelPars$sim.mat))), ncol(ymat))
            aux.list[cases] <- pred.list
            pred.list <- aux.list
            aux.list <- NULL
      }
      x.obs <- getCoordinates(obs)$x
      y.obs <- getCoordinates(obs)$y
      aux.list <- lapply(1:length(modelPars$sim.mat), function (i) {
            aux <- lapply(1:length(pred.list), function(j) {
                  pred.list[[j]][ ,i]
            })
            aux.mat <- do.call("cbind", aux)
            out <- if (!modelPars$stations) {
                  mat2Dto3Darray(aux.mat, x.obs, y.obs)
            } else {
                  aux.mat
            } 
            aux.mat <- NULL
            return(out)
      })
      # Data array - rename dims
      dimNames <- renameDims(obs, modelPars$multi.member)
      obs$Data <- drop(unname(do.call("abind", c(aux.list, along = -1))))
      aux.list <- NULL
      attr(obs$Data, "dimensions") <- dimNames
      attr(obs$Data, "downscaling:method") <- "glm"
      attr(obs$Data, "downscaling:simulation_data") <- modelPars$sim.dataset
      attr(obs$Data, "downscaling:n_eofs") <- n.eofs
      
      # Date replacement
      obs$Dates <- modelPars$sim.dates 
      message("[", Sys.time(), "] Done.")
      return(obs)
}
# End
