#' @title GLM downscaling
#' @description Implementation of GLM's to downscale precipitation data
#' @param y A grid or station data containing the observed climate data for the training period
#' @param modelPars Output object from function ppModelSetup containing the predictors and test data
#' @param simulate Logical. Default is TRUE.
#' @param wet.threshold Value below which precipitation amount is considered zero 
#' @param n.eofs Integer indicating the number of EOFs to be used as predictors 
#' @family downscaling
#' @author J Bedia, M Iturbide
#' @importFrom abind abind
#' @keywords internal

glimpr <- function(y = y, modelPars = modelPars, simulate = TRUE, return.models = FALSE, wet.threshold = wet.threshold, n.pcs = n.pcs) {
      #modelPars <- ppModelSetup(y, x, sim)
      if (is.null(n.pcs)) {
            n.pcs <-  dim(modelPars$pred.mat)[2]
      }
      #pred <- NULL
      #sim <- NULL
      # Prepare predictand matrix
      ymat <- if (isTRUE(modelPars$stations)) {
            y$Data
      } else {
            if (is.null(dim(y$Data))) {
                  as.matrix(y$Data)
            } else {
                  array3Dto2Dmat(y$Data)
            }
      }
      # binarize      
      ymat.bin <- ymat
      ymat.bin[which(ymat < wet.threshold)] <- 0
      ymat.bin[which(ymat >= wet.threshold)] <- 1L
      #obtaining climatological frequencies of occurrences for probability threshold calibration
      wet.prob <- apply(ymat.bin, 2, function(x) {sum(x, na.rm = TRUE) / length(na.exclude(x))})
      # Filter empty gridcell series
      rm.ind <- which(is.na(wet.prob))
      # Valid columns or gridcells (cases)
      cases <- if (length(rm.ind) > 0) {
            (1:ncol(ymat))[-rm.ind]      
      } else {
            1:ncol(ymat)
      }
      # Models
      message("[", Sys.time(), "] Fitting models ...")
      err <- numeric()
      pred.list <- list()
      mod.bin <- list()
      mod.gamma <- list()  
      suppressWarnings(
           for (x in 1:length(cases)){
            mod.bin.x <- glm(ymat.bin[ ,cases[x]] ~ ., data = as.data.frame(modelPars$pred.mat)[,1:n.pcs], family = binomial(link = "logit"))           
            sims.bin <- sapply(1:length(modelPars$sim.mat), function(i) {
                        pred <- predict(mod.bin.x, newdata = as.data.frame(modelPars$sim.mat[[i]])[,1:n.pcs], type = "response")
                        if(simulate == TRUE){
                              rnd <- runif(length(pred), min = 0, max = 1)
                              ind01 <- which(pred > rnd)
                              pred[ind01] <- 1
                              pred[-ind01] <- 0
                        }else{
#                               pred[which(pred >= quantile(pred, 1 - wet.prob[cases[x]], na.rm = TRUE))] <- 1L
                              pred[which(pred >= 0.5)] <- 1L
                              pred <- as.integer(pred)
                        }
                        return(pred)
                  })
            
            if(return.models == FALSE)  mod.bin <- NULL
            wet <- which(ymat[ ,cases[x]] >= wet.threshold)
            # TODO: handle exception when model is over-specified (https://stat.ethz.ch/pipermail/r-help/2000-January/009833.html)
            
            Npcs <- n.pcs+1
            mod.gamma.x = NULL
            
            ###########
            
            while(is.null(mod.gamma.x) & Npcs > 1) {
                  Npcs <-Npcs -1
#                   print(Npcs)
                  mod.gamma.x <- tryCatch(expr = {
                        glm(ymat[ ,cases[x]] ~., data = as.data.frame(modelPars$pred.mat)[ ,1:Npcs], family = Gamma(link = "log"), subset = wet)
                  },
                  error = function(err) {                
                        mod.gamma.x <- NULL
                  })
            }
            
           
            if(length(mod.gamma.x$model)-1 < n.pcs){
                  err[x] <- 1 
            }else{
                  err[x] <- 0
            }
            
            
            ###########
            
            if (is.null(mod.gamma.x)){
                  sims <-  sapply(1:length(modelPars$sim.mat), function(i) {
                        matrix(NA, ncol = 1, nrow = nrow(modelPars$sim.mat[[i]]))
                  })
            }else{
                  sims <- sapply(1:length(modelPars$sim.mat), function(i) {
                        if(simulate == TRUE){
                              predg <- unname(predict(mod.gamma.x, newdata = as.data.frame(modelPars$sim.mat[[i]])[,1:n.pcs], type = "response"))
#                               predg[which(predg<=wet.threshold)] <- 0
                              rgamma(n = length(predg), shape = 1/summary(mod.gamma.x)$dispersion, scale = summary(mod.gamma.x)$dispersion * predg) * sims.bin[,i] 
                        }else{
                               unname(predict(mod.gamma.x, newdata = as.data.frame(modelPars$sim.mat[[i]])[,1:n.pcs], type = "response") ) * sims.bin[,i] 
                        }
                  })
            }
      if(return.models == FALSE)  mod.gamma.x <- NULL
      
      pred.list[[x]] <- unname(sims)
      sims <- NULL
      mod.bin[[x]] <- mod.bin.x 
      mod.gamma[[x]] <- mod.gamma.x
      })
      
      mod <- list("binary" = mod.bin, "gamma" = mod.gamma)
      
      n.err <- length(which(err==1))

      if(n.err > 0){
      n <- round(n.err/length(err)*100, digits = 2)
      warning("In ", n,"% of the locations: glm.fit convergence reached after variable number reduction.")
      }
      

#       err <- lapply(1:length(pred.list), function(u){
#             is.logical(pred.list[[u]])
#       })
#       
#       n.err <- length(which(do.call("abind", err)==TRUE))
#       
#       if(n.err > 0){
#             warning("glm.fit did not converge ", n.err, " times, NA returned.")
#       }
      
      # Refill with NaN series
      #######
      if (length(rm.ind) > 0) {
            aux.list <- rep(list(matrix(nrow = nrow(modelPars$sim.mat[[1]]), ncol = length(modelPars$sim.mat))), ncol(ymat))
            aux.list[cases] <- pred.list
            pred.list <- aux.list
            aux.list <- NULL
      }
      
      x.obs <- getCoordinates(y)$x
      y.obs <- getCoordinates(y)$y
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
      
      #       obs$Data <- 
      o <-drop(unname(do.call("abind", c(aux.list, along = -1))))
      aux.list <- NULL
      # Date replacement
      #       obs$Dates <- modelPars$sim.dates 
      message("[", Sys.time(), "] Done.")
      #       return(obs)
      if(return.models == TRUE) o <- list(o, "models" = mod)
 return(o)
# return(sims.bin)
}
# End
