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
#' @seealso 
#' downscale.train for training a downscaling model
#' prepareNewData for predictor preparation with new (test) data
#' downscale.cv for automatic cross-validation 
#' \href{https://github.com/SantanderMetGroup/downscaleR/wiki/training-downscaling-models}{downscaleR Wiki} for downscaling seasonal forecasting and climate projections.
#' @author J. Bano-Medina
#' @family downscaling.functions
#' @export
#' @examples 
#' # Loading data
#' library(downscaleR)
#' data("VALUE_Iberia_tas")
#' y <- VALUE_Iberia_tas
#' data("NCEP_Iberia_hus850", "NCEP_Iberia_psl", "NCEP_Iberia_ta850")
#' x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_psl, NCEP_Iberia_ta850)
#' 
#' # Example1: Basic example (without cross-validation)
#' data <- prepareData(x = x, y = y, spatial.predictors = list(v.exp = 0.95))
#' model.analogs <- downscale.train(data, method = "analogs", n.analogs = 1)
#' newdata <- prepareNewData(x,data)
#' pred <- downscale.predict(newdata, model.analogs)
#' # This produces the same result as model.analogs$pred
#' 
#' # Example2:  Splitting data in train and test (simple cross-validation)
#' xT <- subsetGrid(x, years = 1983:1999)  # training predictors
#' yT <- subsetGrid(y, years = 1983:1999)   # training predictands
#' data <- prepareData(xT,yT)       # preparing the data
#' model.analogs <- downscale.train(data, method = "analogs", n.analogs = 1)
#' xt <- subsetGrid(x, years = 2000)       # test predictors
#' yt <- subsetGrid(y, years = 2000)       # test predictors
#' newdata <- prepareNewData(xt,data)     # preparing the new predictors
#' pred  <- downscale.predict(newdata, model.analogs)  # predicting
#' # Plotting the results for station 5
#' plot(yt$Data[,5],pred$Data[,5])

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
        
        if (is.null(model$model$atomic_model[[i]])) {
          yp[,i] <- rep(NaN,n.obs)  
        }
        else {
          if (model$model$method == "analogs") {model$model$atomic_model[[i]]$dates$test <- getRefDates(newdata)}
          if (is.null(model$model$atomic_model[[i]])) {
            yp[,i] = rep(NA, 1, n.obs)
          } else {
            yp[,i] <- downs.predict(xx, model$model$method, model$model$atomic_model[[i]])
          }
        }
      }
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
        if (is.null(model$model$atomic_model[[i]])) {
          yp[,i] <- rep(NaN,n.obs)  
        }
        else {
          if (model$model$method == "analogs") {model$model$atomic_model[[i]]$dates$test <- getRefDates(newdata)}
          if (is.null(model$model$atomic_model[[i]])) {
            yp[,i] = rep(NA, 1, n.obs)
          } else {
            yp[,i] <- downs.predict(xx, model$model$method, model$model$atomic_model[[i]])
          }
        }
      }
    }
    if (isRegular(model$pred)) {
      yp <- mat2Dto3Darray(yp, x = model$pred$xyCoords$x, y = model$pred$xyCoords$y)
    }
    return(yp)
  })
  
  if (isRegular(model$pred)) {
    if (n > 1) {
      p <- array(unlist(p), dim = c(dim(p[[1]]),n)) %>% aperm(c(4,1,2,3))
      dimNames <- getDim(redim(model$pred))}
    else {
      p <- array(unlist(p), dim = dim(p[[1]]))
      dimNames <- getDim(model$pred)}
  }
  else {
    if (n > 1) {
      p <- array(unlist(p), dim = c(dim(p[[1]]),n)) %>% aperm(c(3,1,2))
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
#' @importFrom deepnet nn.predict
#' @author J. Bano-Medina
downs.predict <- function(x, method, atomic_model){
  switch(method,
         analogs = pred <- analogs.test(x, atomic_model$dataset_x, atomic_model$dataset_y, atomic_model$dates, atomic_model$info),
         GLM     = pred <- glm.predict(x, atomic_model$weights, atomic_model$info),
         NN      = pred <- nn.predict(atomic_model, x)) 
  return(pred)}