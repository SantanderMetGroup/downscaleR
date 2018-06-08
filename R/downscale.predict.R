##############################################################################################################
#                     GENERAL DOWNSCALING                                                                    #
##############################################################################################################
##     downscale.predict.R Downscale climate data for a given statistical model.
##
##     Copyright (C) 2017 Santander Meteorology Group (http://www.meteo.unican.es)
##
##     This program is free software: you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 3 of the License, or
##     (at your option) any later version.
## 
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
## 
##     You should have received a copy of the GNU General Public License
##     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Downscale climate data for a given statistical model.
#' @description Downscale data to local scales by statistical models previously obtained by \code{\link[downscaleR]{downscale.train}}.
#' @param newdata The grid data. It should be an object as returned by  \code{\link[downscaleR]{prepareNewData}}.
#' @param model An object containing the statistical model as returned from  \code{\link[downscaleR]{downscale.train}}.
#' @return A regular/irregular grid object.
#' @author J. Bano-Medina
#' @export
#' @examples 
#' # Loading predictors
#' x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_ta850)
#' x <- subsetGrid(x, years = 1985:1995)
#' # Loading predictands
#' y <- VALUE_Iberia_pr
#' y <- getTemporalIntersection(obs = y,prd = x, "obs" )
#' x <- getTemporalIntersection(obs = y,prd = x, "prd" )
#' ybin <- binaryGrid(y, threshold = 1)
#' x <- scaleGrid(x, type = "standardize")
#' # Prepare predictors and predictands
#' xyT     <- prepareData(x = x, y = y)
#' xyT.bin <- prepareData(x = x, y = ybin)
#' xyt     <- prepareNewData(newdata = x, data.structure = xyT)
#' xyt.bin <- prepareNewData(newdata = x, data.structure = xyT.bin)
#' # Downscaling PRECIPITATION
#' # ... via analogs ...
#' model <- downscale.train(xyT, method = "analogs",
#'                          sel.fun = "mean")
#' pred <- downscale.predict(xyt, model = model)
#' # ... via a logistic regression (ocurrence of precipitation)
#' # and gaussian regression (amount of precipitation) ...
#' model.ocu <- downscale.train(xyT.bin, method = "GLM",
#'                              family = binomial(link = "logit"))
#' pred <- downscale.predict(xyt, model = model.ocu)
#' model.reg <- downscale.train(xyT, method = "GLM",
#'                              family = "gaussian", filter = ">0")
#' pred <- downscale.predict(xyt, model = model.reg)
#' # ... via a neural network ...
#' model.ocu <- downscale.train(xyT.bin, method = "NN",
#'                              learningrate = 0.1, numepochs = 10, hidden = 5,
#'                              output = 'linear')
#' pred <- downscale.predict(xyt, model = model.ocu)
#' model.reg <- downscale.train(xyT, method = "NN",
#'                              learningrate = 0.1, numepochs = 10,
#'                              hidden = 5, output = 'linear')
#' pred <- downscale.predict(xyt, model = model.reg)
#' # Downscaling PRECIPITATION - Local model with the closest
#' # 4 grid points and multisite linear regression.
#' xyT.local <- prepareData(x = x, y = y,
#'                          local.predictors = list(vars = "hus@850",n = 4))
#' xyt.local <- prepareNewData(newdata = x, data.structure = xyT.local)
#' model <- downscale.train(xyT.local,method = "analogs")
#' pred <- downscale.predict(xyt.local, model = model)
#' # Downscaling PRECIPITATION - Principal Components (PCs)
#' # and gamma regression for the amount of precipitation
#' xyT.pc <- prepareData(x = x,y = y,
#'                       spatial.predictors = list(which.combine = getVarNames(x),v.exp = 0.9))
#' xyt.pc <- prepareNewData(newdata = x, data.structure = xyT.pc)
#' model <- downscale.train(xyT.pc, method = "analogs")
#' pred <- downscale.predict(xyt.pc, model = model)

downscale.predict <- function(newdata, model) {
  n <- length(newdata$x.global) # number of members
  p <- lapply(1:n, function(z) {
    # Multi-site
    if (model$model$site == "multi") {
      xx <- newdata$x.global[[z]]
      if (model$model$method == "analogs") {model$model$atomic_model$dates$test <- getRefDates(newdata)}
      yp <- as.matrix(downs.predict(xx, model$model$method, model$model$atomic_model))}
    # Single-site
    else if (model$model$site == "single") {
      stations <- length(model$model$atomic_model)
      if (attr(newdata$x.global,"nature") == "local") {
        n.obs <- nrow(newdata$x.local[[1]][[z]])}
      else {
        n.obs <- nrow(newdata$x.global[[z]])}
      yp <- array(data = NA, dim = c(n.obs,stations))
      for (i in 1:stations) {
        if (!is.null(newdata$x.local)) {
          xx = newdata$x.local[[i]][[z]]}
        else {
          xx <- newdata$x.global[[z]]}
        if (model$model$method == "analogs") {model$model$atomic_model[[i]]$dates$test <- getRefDates(newdata)}
        yp[,i] <- downs.predict(xx, model$model$method, model$model$atomic_model[[i]])}
    }
    # Mix - Global predictors with local predictors
    else if (model$model$site == "mix") {
      stations <- length(model$model$atomic_model)
      if (attr(newdata$x.global,"nature") == "local") {
        n.obs <- nrow(newdata$x.local[[1]][[z]])}
      else {
        n.obs <- nrow(newdata$x.global[[z]])}
      yp <- array(data = NA, dim = c(n.obs,stations))
      for (i in 1:stations) {
        xx1 = newdata$x.local[[i]][[z]]
        xx2 = newdata$x.global[[z]]
        xx <- cbind(xx1,xx2)
        if (model$model$method == "analogs") {model$model$atomic_model[[i]]$dates$test <- getRefDates(newdata)}
        yp[,i] <- downs.predict(xx, model$model$method, model$model$atomic_model[[i]])}
    }
    if (isRegular(model$pred)) {
      yp <- mat2Dto3Darray(yp, x = model$pred$xyCoords$x, y = model$pred$xyCoords$y)
    }
    return(yp)
  })
  
  if (isRegular(model$pred)) {
    if (n > 1) {
      p <- array(unlist(p), dim = c(n,dim(p[[1]])))
      dimNames <- getDim(redim(model$pred))}
    else {
      p <- array(unlist(p), dim = dim(p[[1]]))
      dimNames <- getDim(model$pred)}
  }
  else {
    if (n > 1) {
      p <- array(unlist(p), dim = c(n,dim(p[[1]])))
      dimNames <- getDim(redim(model$pred,loc = TRUE))}
    else {
      p <- array(unlist(p), dim = dim(p[[1]]))
      dimNames <- getDim(model$pred)}
  }
  pred <- model$pred
  pred$Data <- p
  attr(pred$Data, "dimensions") <- dimNames
  pred$Dates <- newdata$Dates
  return(pred)
}

##############################################################################################################
#                     DOWNSCALING                                                                            #
##############################################################################################################
#' @title Switch to selected downscale method.
#' @description Internal function of \code{\link[downscaleR]{downscale.predict}} that switches to the corresponding method.
#' @param x The grid data. Class: matrix.
#' @param method The method of the given model.
#' @param atomic_model An object containing the statistical model of the selected method.
#' @return A matrix with the predictions.
#' @details This function is internal and should not be used by the user. The user should use \code{\link[downscaleR]{downscale.predict}}.
#' @author J. Bano-Medina
#' @export
downs.predict <- function(x, method, atomic_model){
  switch(method,
         analogs = pred <- analogs.test(x, atomic_model$dataset_x, atomic_model$dataset_y, atomic_model$dates, atomic_model$info),
         GLM     = pred <- glm.predict(x, atomic_model$weights, atomic_model$info),
         NN      = pred <- nn.predict(atomic_model, x)) 
  return(pred)}