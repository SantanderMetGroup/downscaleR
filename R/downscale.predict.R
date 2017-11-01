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
#' @param newdata The grid data. It should be an object as returned by  \code{\link[downscaleR]{prepare_newdata}}.
#' @param model An object containing the statistical model as returned from  \code{\link[downscaleR]{downscale.train}}.
#' @return An object with the predictions.
#' @details The function can downscale in both global and local mode, though not simultaneously.
#' @author J. Bano-Medina
#' @export
#' @importFrom MASS ginv 
#' @importFrom matlab reshape repmat
#' @examples 
#' # Loading predictors
#' x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_ta850)
#' x <- subsetGrid(x, years = 1985:1995)
#' # Loading predictands
#' y <- VALUE_Iberia_pr
#' y <- getTemporalIntersection(obs = y,prd = x, "obs" )
#' ybin <- convert2bin(y, threshold = 1)
#' # Prepare predictors
#' xT <- prepare_predictors(x = x, y = y)
#' # Downscaling PRECIPITATION
#' #' # ... via analogs ...
#' model.ocu <- downscale.train(ybin, xT, method = "analogs", sel.fun = "mean")
#' model.reg <- downscale.train(y, xT, method = "analogs", sel.fun = "mean")
#' pred.ocu <- downscale.predict(xT, model.ocu)
#' pred.reg <- downscale.predict(xT, model.reg)
#' # ... via a linear model ...
#' model.ocu <- downscale.train(ybin, xT, method = "GLM" ,family = binomial(link = "logit"))
#' model.reg <- downscale.train(y, xT, method = "MLR", fitting = "MP")
#' pred.ocu <- downscale.predict(xT,model.ocu)
#' pred.reg <- downscale.predict(xT,model.reg)
#' pred <- pred.ocu$Data * pred.reg$Data
#' # ... via a extreme learning machine ...
#' model.ocu <- downscale.train(ybin, xT, method = "ELM", neurons = 200)
#' model.reg <- downscale.train(y, xT, method = "ELM", neurons = 200)
#' pred.ocu <- downscale.predict(xT,model.ocu)
#' pred.reg <- downscale.predict(xT,model.reg)
#' pred <- pred.ocu$Data * pred.reg$Data
#' # ... via a neural network ...
#' model.ocu <- downscale.train(ybin, xT, method = "NN", learningrate = 0.1, numepochs = 10, hidden = 5, output = 'linear') 
#' model.reg <- downscale.train(y, xT, method = "NN", learningrate = 0.1, numepochs = 10, hidden = 5, output = 'linear') 
#' pred.ocu <- downscale.predict(xT,model.ocu)
#' pred.reg <- downscale.predict(xT,model.reg)
#' pred <- pred.ocu$Data * pred.reg$Data
#' # Downscaling PRECIPITATION - Local model with the closest 4 grid points.
#' xT <- prepare_predictors(x = x,y = y,local.predictors = list(neigh.vars = "shum850",n.neighs = 4))
#' model.ocu <- downscale.train(ybin, xT, method = "MLR", fitting = 'MP')
#' model.reg <- downscale.train(y, xT, method = "MLR", fitting = 'MP')
#' pred.ocu <- downscale.predict(xT,model.ocu)
#' pred.reg <- downscale.predict(xT,model.reg)
#' pred <- pred.ocu$Data * pred.reg$Data
#' # Downscaling PRECIPITATION - Principal Components (PCs)
#' xT <- prepare_predictors(x = x,y = y, PCA = list(which.combine = getVarNames(x),v.exp = 0.9))
#' model.ocu <- downscale.train(ybin, xT, method = "MLR" ,fitting = 'MP')
#' model.reg <- downscale.train(y, xT, method = "MLR" ,fitting = 'MP')
downscale.predict <- function(newdata, model) {
  dimNames <- getDim(model$pred)
  pred <- model$pred
  pred$Dates <- newdata$Dates
  if (!isTRUE(model$conf$singlesite)) {
    if (!is.list(newdata$x.global)) {
      xx <- newdata$x.global}
    else {
      xx <- newdata$x.global$member_1}
    pred$Data <- downs.predict(xx, model$conf$method, model$conf$atomic_model)}
  else {
    if (is.null(dim(model$pred$Data))) {stations <- 1; n.obs <- length(model$pred$Data)}
    else{stations <- dim(model$pred$Data)[2]; n.obs <- length(newdata$Dates$start)}
    if (is.null(newdata$Dates)) {n.obs <- length(newdata$y$Dates$start)}
    pred$Data <- array(data = NA, dim = c(n.obs,stations))
    for (i in 1:stations) {
      if (!is.null(newdata$x.local)) {
        if (!is.list(newdata$x.local[[i]])) {
          xx = newdata$x.local[[i]]}
        else {
          xx = newdata$x.local[[i]]$member_1}}
      else {
        if (!is.list(newdata$x.global)) {
          xx <- newdata$x.global}
        else {
          xx <- newdata$x.global$member_1}}
      pred$Data[,i] <- downs.predict(xx, model$conf$method, model$conf$atomic_model[[i]])}}
  attr(pred$Data, "dimensions") <- dimNames
  return(pred)}

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
  if (!is.null(atomic_model$info$fitting) && (atomic_model$info$fitting == "MP" || atomic_model$info$fitting == "MP+L2")) {
    x <- cbind(x,array(data = 1,dim = c(nrow(x), 1)))}
  switch(method,
         analogs = pred <- analogs.test(x, atomic_model$dataset_x, atomic_model$dataset_y, atomic_model$info),
         GLM     = pred <- glm.predict(x, atomic_model),
         MLR     = pred <- mlr.predict(x, atomic_model$weights, atomic_model$info),
         ELM     = pred <- elm.predict(x, atomic_model$weights, atomic_model$info),
         NN      = pred <- nn.predict(atomic_model, x)) 
  return(pred)}