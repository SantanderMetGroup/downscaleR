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
#' x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_psl, NCEP_Iberia_ta850)
#' x <- subsetGrid(x, years = 1985:1995)
#' # Loading predictands
#' y <- VALUE_Iberia_pr
#' # Prepare predictors and prepare new data
#' xtrain <- prepare_predictors(x = x, y = y)
#' xtest <- prepare_newdata(newdata = x, predictor = xtrain)
#' # Downscaling a generalized linear model. In particular a logistic regression.
#' # model <- downscale.train(ybin, xtrain, method = "GLM" ,family = binomial(link = "logit"))
#' # str(model)
#' # pred <- downscale.predict(xtest, model)
#' # str(pred)
#' # Downscaling a linear model..
#' model <- downscale.train(xtrain$y, xtrain, method = "MLR", fitting = 'MP')
#' str(model)
#' pred <- downscale.predict(xtest, model)
#' str(pred)
#' # Downscaling a non-linear model via a regularized L2 Extreme learning Machine..
#' model <- downscale.train(xtrain$y, xtrain, method = "ELM", fitting = 'MP+L2')
#' str(model)
#' pred <- downscale.predict(xtest, model)
#' str(pred)
#' # Downscaling a non-linear model via Neural Networks..
#' model <- downscale.train(xtrain$y, xtrain, method = "NN", learningrate = 0.1, numepochs = 10, hidden = 5, output = 'linear') 
#' str(model)
#' pred <- downscale.predict(xtest, model)
#' str(pred)
#' # Downscaling a local model with the closest 4 grid points.
#' xtrain <- prepare_predictors(x = x,y = y,local.predictors = list(neigh.vars = "shum850",n.neighs = 4))
#' xtest <- prepare_newdata(newdata = x, predictor = xtrain)
#' model <- downscale.train(xtrain$y, xtrain, method = "MLR", fitting = 'MP')
#' str(model)
#' pred <- downscale.predict(xtest, model)
#' str(pred)
#' # Donwscaling with Principal Components (PCs)
#' xtrain <- prepare_predictors(x = x,y = y, PCA = list(v.exp = 0.9))
#' xtest <- prepare_newdata(newdata = x, predictor = xtrain)
#' model <- downscale.train(xtrain$y, xtrain, method = "ELM" ,fitting = 'MP')
#' str(model)
#' pred <- downscale.predict(xtest, model)
#' str(pred)
downscale.predict <- function(newdata, model) {
  dimNames <- getDim(model$predictors$y)
  pred <- model$pred
  pred$Dates <- newdata$Dates
  if (is.null(newdata$newdata.local)) {
    pred$Data <- downs.predict(newdata$newdata.global$member_1, model$conf$method, model$conf$atomic_model)}
  if (!is.null(newdata$newdata.local)) {
    stations  <- length(newdata$newdata.local)
    pred$Data <- array(data = NA, dim = c(nrow(newdata$newdata.local[[1]]$member_1), stations))
    for (i in 1:stations) {
      pred$Data[,i] <- downs.predict(newdata$newdata.local[[i]]$member_1, model$conf$method, model$conf$atomic_model[[i]])}}
  attr(pred$Data, "dimensions") <- dimNames
  pred <- list("pred" = pred)
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
         GLM   = pred <- glm.predict(x, atomic_model),
         MLR   = pred <- mlr.predict(x, atomic_model$weights, atomic_model$info),
         ELM   = pred <- elm.predict(x, atomic_model$weights, atomic_model$info),
         NN    = pred <- nn.predict(atomic_model, x)) 
  return(pred)}

##############################################################################################################
#                     Generalized Linear Models (GLM)                                                        #
##############################################################################################################
#' @title Donwscaling with a given generalized linear model (GLM).
#' @description Donwscaling with a given generalized linear models (GLM) calculated in \code{\link[downscaleR]{downscale.predict}} via \code{\link[downscaleR]{glm.train}}.
#' @param x The grid data. Class: matrix.
#' @param weights Object as returned from \code{\link[downscaleR]{glm.train}}
#' @return The predicted matrix.
#' @details Predicts by using the base function \code{\link[stats]{predict}}. This function is internal and should not be used by the user. The user should use \code{\link[downscaleR]{downscale.predict}}.
#' @author J. Bano-Medina
#' @export
glm.predict <- function(x, weights){
  pred <- array(data = NA , dim = c(dim(x)[1],length(weights)))
  for (i in 1:length(weights)) {
    df <- data.frame(x); colnames(df) <- paste0('X',2:106)
    if (weights[[i]]$family$family == "binomial") {
      pred[,i] <- predict(weights[[i]], newdata = df, type = 'response')}
    else{
      pred[,i] <- predict(weights[[i]], newdata = df)}}
  return(pred)}

##############################################################################################################
#                     Multiple Linear Regression (MLR)                                                       #
##############################################################################################################
#' @title Donwscaling with a given multiple linear regression model (MLR).
#' @description Donwscaling with a given multiple linear regression (MLR) calculated in \code{\link[downscaleR]{downscale.predict}} via \code{\link[downscaleR]{mlr.train}}.
#' @param x The grid data. Class: matrix.
#' @param weights Object as returned from \code{\link[downscaleR]{mlr.train}}
#' @param info  This is a list containing the information/parameters used in \code{\link[downscaleR]{mlr.train}}. Basically is used to distinguish between Least Squares (LS) and Moore-Penrose (MP).
#' @return The predicted matrix.
#' @details Predicts by using the base function \code{\link[stats]{predict}}. This function is internal and should not be used by the user. The user should use \code{\link[downscaleR]{downscale.predict}}.
#' @author J. Bano-Medina
#' @export
mlr.predict <- function(x, weights, info){
  if (info$fitting == 'LS') {
    pred <- array(data = NA , dim = c(dim(x)[1],length(weights)))
    for (i in 1:length(weights)) {
      df <- data.frame(x); colnames(df) <- paste0('X',2:106)
      pred[,i] <- predict(weights[[i]], newdata = df, type = 'response')}}
  if (info$fitting == 'MP' || info$fitting == 'MP+L2') {
    pred <- data.matrix(x) %*% weights}
  return(pred)}

##############################################################################################################
#                     Extreme Learning Machine (ELM)                                                         #
##############################################################################################################
#' @title Donwscaling with a given extreme learning machine model (ELM).
#' @description Donwscaling with a given extreme learning machine (ELM) calculated in \code{\link[downscaleR]{downscale.predict}} via \code{\link[downscaleR]{elm.train}}.
#' @param x The grid data. Class: matrix.
#' @param weights Matrix of coefficients as returned from \code{\link[downscaleR]{elm.train}}
#' @param info  This is a list containing the information/parameters used in \code{\link[downscaleR]{elm.train}}. Containes the optional parameters: fitting, neurons, Act.F, area.region and area.module. See \code{\link[downscaleR]{elm.train}} for information relevant to these parameters.
#' @return The predicted matrix.
#' @details This function is internal and should not be used by the user. The user should use \code{\link[downscaleR]{downscale.predict}}.
#' @author J. Bano-Medina
#' @export
elm.predict <- function(x, weights, info){
  h <- activation(data.matrix(x) %*% weights[[1]], info$Act.F)
  pred <- h %*% weights[[2]]
  return(pred)}