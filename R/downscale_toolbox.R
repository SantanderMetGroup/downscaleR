##############################################################################################################
#                     GENERAL DOWNSCALING                                                                    #
##############################################################################################################
##     downscale.train.R Downscale climate data.
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

#' @title Downscale climate data.
#' @description Downscale data to local scales by statistical methods: generalized linear models (GLM), multiple linear regression (MLR), Extreme Learning Machine (ELM) and Neural Networks (NN). 
#' @param predictands The observations dataset. It should be an object as returned by \pkg{loadeR}.
#' @param predictors The input grid as returned by \code{\link{prepare_predictors}}.
#' @param method Type of transer function. Options are: GLM, MLR, ELM and NN. The default is MLR.
#' @param ... Optional parameters required for the specific method chosen. 
#' GLM: Uses the glm base function. Parameters: family. 
#' MLR: Uses wether an internal function or the lm base function. Parameters: fitting = c("LS","MP","MP+L2"). You can fit the model via the Least Squares (LR) lm R function, or via the Moore-Penrose (MP) inverse, or via the Moore-Penrose regularized by a L2 penalty (MP+L2)
#' ELM: Uses an internal function. Parameters: fitting = c("MP","MP+L2"), neurons = 100, Act.F = 'sig', area.region = NULL, area.module = NULL. 
#' The number of hidden neurons is selected in the parameter "neurons" along with its activation function in "Act.F.
#' For a modular ELM, called Receptive Fields ELM (RF-ELM), special parameters are needed: area.region and area.module. 
#' The area.region parameter means the area covered by the grid (i.e, for the following gridpoints: (46,2), (48,2), (50,2), (46,4), (48,4), (50,4); your area is c(3,2)). 
#' The area.module parameter means the size of the square mask (i.e., if area.module = 2, this means that a hidden neuron will be fed by all variables included in a random selected square of 2 consecutive gridpoints in latitude and 2 consecutive points in longitude)  
#' NN:  Uses the nn.train function of the deepnet library. Parameters: initW = NULL, initB = NULL, hidden = c(10), activationfun = "sigm", learningrate = 0.8, momentum = 0.5, learningrate_scale = 1, output = "sigm", numepochs = 3, batchsize = 100, hidden_dropout = 0, visible_dropout = 0
#' @return A list of objects that contains the model, the prediction on the train dataset, the predictors and predictands used.
#' @details The function can downscale in both global and local mode, though not simultaneously. 
#' If there is perfect collinearity among predictors, then the matrix will not be invertible and the downscaling will fail.
#' 
#' 
#' Appends the attribute \code{subset = "convert2bin"} in the \code{Variable} element.
#' 
#' @author J. Bano-Medina
#' @export
#' @importFrom MASS ginv 
#' @family subsetting
#' @examples
#' data("VALUE_Iberia_tp")
#' # Take a look at the data:
#' head(VALUE_Iberia_tp$Data)
#' # Convert to complete binary variable:
#' bin.total <- convert2bin(VALUE_Iberia_tp,threshold = 1)
#' head(bin.total$Data)
#' # Convert to partial binary variable:
#' bin.partial <- convert2bin(VALUE_Iberia_tp,threshold = 1, partial = TRUE)
#' head(bin.partial$Data)
##############################################################################################################
#                     GENERAL DOWNSCALING                                                                    #
##############################################################################################################
downscale.train <- function(predictands, predictors, method, ...) {
  dimNames <- getDim(predictors$y)
  pred <- predictors$y
  if (is.null(predictors$x.local)) {
    atomic_model <- downs.train(predictors$x.global, predictands$Data, method, ...)
    pred$Data    <- downs.predict(predictors$x.global, method, atomic_model)}
  if (!is.null(predictors$x.local)) {
    stations      <- length(predictors$x.local)
    atomic_model  <- vector("list",stations)
    pred$Data     <- array(data = NA, dim = c(nrow(predictors$x.local[[1]]$member_1),stations))
    for (i in 1:stations) {
      atomic_model[[i]] <- downs.train(predictors$x.local[[i]]$member_1, predictands$Data[,i], method, ...)
      pred$Data[,i]     <- downs.predict(predictors$x.local[[i]]$member_1, method, atomic_model[[i]])}}
  attr(pred$Data, "dimensions") <- dimNames
  model <- list("pred" = pred, "conf" = list("method" = method, "atomic_model" = atomic_model), "predictors" = predictors, "predictands" = predictands)
  return(model)}

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
downs.train <- function(x, y, method = "MLR", ...) {
  switch(method,
         GLM   = atomic_model <- glm.train(x, y, ...),
         MLR   = atomic_model <- mlr.train(x, y, ...),
         ELM   = atomic_model <- elm.train(x, y, ...),
         NN    = atomic_model <- nn.train(x, y,  ...))
  return(atomic_model)}

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
glm.train <- function(x, y, family) {
  neurons <- vector("integer",2)
  neurons[1] <- dim(x)[2] + 1
  if (is.null(dim(y)[2])) {neurons[2] <- 1}
  else{neurons[2] <- dim(y)[2]}
  weights <- vector("list",neurons[2])
  for (i in 1:neurons[2]) {
    df <- data.frame(cbind(y[,i],x))
    weights[[i]] <- glm(X1~.,data = df, family)}
  return(weights)}

glm.predict <- function(x, weights){
  pred <- array(data = NA , dim = c(dim(x)[1],length(weights)))
  for (i in 1:length(weights)) {
    df <- data.frame(x); colnames(df) <- paste0('X',2:106)
    pred[,i] <- predict.glm(weights[[i]],df)}
  return(pred)}
##############################################################################################################
#                     Multiple Linear Regression (GLM)                                                       #
##############################################################################################################





mlr.train <- function(x, y, fitting = 'MP'){
  neurons <- vector("integer",2)
  neurons[1] <- dim(x)[2] + 1
  if (is.null(dim(y)[2])) {neurons[2] <- 1}
  else{neurons[2] <- dim(y)[2]}
  if (fitting == 'LS') {
    weights <- vector("list",neurons[2])
    for (i in 1:neurons[2]) {
      df <- data.frame(cbind(y[,i],x))
      weights[[i]] <- lm(X1~.,data = df)}}
  if (fitting == 'MP') {
    weights <- array(data = NaN, c(neurons[1], neurons[2]))
    x <- cbind(x,array(data = 1,dim = c(nrow(x), 1)))
    weights <- ginv(data.matrix(x)) %*% data.matrix(y)}
  if (fitting == 'MP+L2') {
    weights <- array(data = NaN, c(neurons[1], neurons[2]))
    x <- cbind(x,array(data = 1,dim = c(nrow(x), 1)))
    C <- c(0.0001,0.001,0.01,0.1,1,10,100,1000)
    error <- vapply(1:length(C), FUN.VALUE = numeric(length = 1),FUN = function(z) {
      cross.v(x, y, "MP+L2",C[z])})
    C_winner <- C[which(error == min(error))]
    print(C_winner)
    weights <- solve((t(data.matrix(x)) %*% data.matrix(x) + (diag(dim(data.matrix(x))[2]) / C_winner))) %*% (t(data.matrix(x)) %*% data.matrix(y))}
  return(list("weights" = weights, "info" = list("fitting" = fitting)))}

mlr.predict <- function(x, weights, info){
  if (info$fitting == 'LS') {
    pred <- array(data = NA , dim = c(dim(x)[1],length(weights)))
    for (i in 1:length(weights)) {
      df <- data.frame(x); colnames(df) <- paste0('X',2:106)
      pred[,i] <- predict.lm(weights[[i]],df)}}
  if (info$fitting == 'MP' || info$fitting == 'MP+L2') {
    pred <- data.matrix(x) %*% weights}
  return(pred)}
##############################################################################################################
#                     Extreme Learning Machine (ELM)                                                         #
##############################################################################################################
elm.train <- function(x, y, fitting = 'MP', neurons = 100, Act.F = 'sig', area.region = NULL, area.module = NULL){
  neurons[2] <- neurons
  neurons[1] <- dim(x)[2] + 1
  if (is.null(dim(y)[2])) {neurons[3] <- 1}
  else{neurons[3] <- dim(y)[2]}
  x   <- cbind(x,array(data = 1,dim = c(nrow(x), 1)))
  weights <- vector('list',2)
  weights[[1]] <- initialize.weights(x, neurons, area.region, area.module)
  h <- activation(x %*% weights[[1]], Act.F)
  if (fitting == 'MP') {
    weights[[2]] <- ginv(h) %*% y}
  if (fitting == 'MP+L2') {
    C <- c(0.0001,0.001,0.01,0.1,1,10,100,1000)
    error <- vapply(1:length(C), FUN.VALUE = numeric(length = 1),FUN = function(z) {
      cross.v(h, y, "MP+L2",C[z])})
    C_winner <- C[which(error == min(error))] ; print(C_winner)
    weights[[2]] <- solve(t(h) %*% h + (diag(dim(h)[2]) / C_winner)) %*% (t(h) %*% data.matrix(y))}
  return(list("weights" = weights, "info" = list("fitting" = fitting, "neurons" = neurons, "Act.F" = Act.F, "area.region" = area.region, "area.module" = area.module)))}

elm.predict <- function(x, weights, info){
  h <- activation(data.matrix(x) %*% weights[[1]], info$Act.F)
  pred <- h %*% weights[[2]]
  return(pred)}

##############################################################################################################
#                     AUXILIARY FUNCTIONS                                                                    #
##############################################################################################################
initialize.weights <- function(x, neurons, area.region = NULL, area.module = NULL) {
  if (is.null(area.region)) {
    #inp.w <- array(data = runif(neurons[1]*neurons[2], min = -1, max = 1), c(neurons[1],neurons[2]))
    #inp.w <- apply(X = inp.w, MARGIN = 2, function(x) {x <- x / norm(x = x, type = "2")})
    r <- 2*((neurons[1]) ** (-0.5))
    inp.w <- array(data = runif(neurons[1]*neurons[2], min = -r, max = r), c(neurons[1],neurons[2]))
    inp.w[neurons[1],] <- runif(neurons[2], min = -2, max = 2)
    return(inp.w)}
  if (!is.null(area.region)) {
    r <- 2*((neurons[1]) ** (-0.5))
    inp.w <- array(data = runif(neurons[1]*neurons[2], min = -r, max = r), c(neurons[2],neurons[1]))
    inp.w[,neurons[1]] <- runif(neurons[2], min = -2, max = 2)
    #inp.w <- array(data = runif(neurons[1]*neurons[2], min = -1, max = 1), c(neurons[2],neurons[1]))
    r1 <- sample(1:(area.region[1] - area.module + 1), neurons[2], replace = T)
    r2 <- sample(1:(area.region[2] - area.module + 1), neurons[2], replace = T)
    mask <- array(data = 0, c(area.region,neurons[2]))
    for (i in 1:neurons[2]) {mask[r1[i]:(r1[i] + area.module - 1),r2[i]:(r2[i] + area.module - 1),i] <- 1}  
    dim(mask) <- c(area.region[1]*area.region[2],neurons[2])
    num_var <- (neurons[1] - 1)/dim(mask)[1]
    inp.w[,1:(neurons[1] - 1)] <- inp.w[,(1:neurons[1] - 1)] * repmat(t(mask),c(1,num_var)); inp.w <- t(inp.w)
    #inp.w <- apply(X = inp.w, MARGIN = 2, function(x) {x <- x / norm(x = x, type = "2")})
    return(inp.w)}}

activation <- function(temp.h, Act.F){
  switch(Act.F,
         lin = temp.h,
         sig = 1/(1 + exp(-temp.h)))}

cross.v <- function(x, y, fitting, C) {
  data  <- sample(x = 1:(floor(nrow(x)/4) * 4),replace = FALSE)
  data  <- reshape(data.matrix(data[1:(floor(nrow(x)/4) * 4)]),c(floor(nrow(x)/4),4))
  error <- lapply(1:4,FUN = function(z) {
    data.train <- data.matrix(x)[data[,-z],]; target.train <- data.matrix(y)[data[,-z],]
    data.test  <- data.matrix(x)[data[, z],]; target.test  <- data.matrix(y)[data[, z],]
    weights <- solve(t(data.train) %*% data.train + (diag(dim(data.train)[2]) / C)) %*% (t(data.train) %*% target.train)
    pred <- mlr.predict(data.test, weights, list(fitting = "MP+L2"))
    error <- sum(apply((pred - target.test) ** 2, MARGIN = 1, FUN = sum))/dim(x)[1]
    return(error)})
  error <- mean(as.numeric(error))
  return(error)}
