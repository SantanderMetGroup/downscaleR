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
#' @param predictors The input grid as returned by \code{\link[downscaleR]{prepare_predictors}}.
#' @param method Type of transer function. Options are: GLM, MLR, ELM and NN. The default is MLR.
#' @param ... Optional parameters. These parameters are different depending on the method selected. Every parameter has a default value set in the atomic functions in case that no selection is wanted. Everything concerning these parameters is explained in the section \code{Details}. However, if wanted, the atomic functions can be seen here: \code{\link[downscaleR]{glm.train}}, \code{\link[downscaleR]{mlr.train}}, \code{\link[downscaleR]{elm.train}} and \code{\link[deepnet]{nn.train}}.  

#' @details The function can downscale in both global and local mode, though not simultaneously. 
#' If there is perfect collinearity among predictors, then the matrix will not be invertible and the downscaling will fail.
#' 
#' \strong{Generalized Linear Models (GLM)}
#' 
#' This function uses \code{\link[stats]{glm}}. The unique optional parameter is \code{family} with default \code{gaussian}. The possible family options are: gaussian, binomial, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, quasipoisson. Indeed, family is a function itself of the form
#' family(object,...) where you can specified a link (i.e., a specification for the model link function), see \code{\link[stats]{family}}. For example for a logistic regression the optional parameter would be family = binomial(link = "logit"). The optional parameters of each method are axplained here:
#' 
#' \strong{Multiple Linear Regression (MLR)}
#' 
#' If you want to downscale by multiple linear regression there is only one optional parameter called, \code{fitting}, with options: fitting = c("LS","MP","MP+L2"). This options refers to Least Squares (LS), Moore-Penrose (MP) and L2 penalty Moore-Penrose (MP+L2). LS uses the \code{\link[stats]{lm}} R function, whereas MP and MP+L2 downscales via an internal function \code{\link[downscaleR]{mlr.train}}. "MP" is the default option.
#' 
#' \strong{Extreme Learning Machine (ELM)}
#' 
#' If you want to downscale via an Extreme Learning Machine there are 5 optional parameters: \code{fitting}, \code{neurons}, \code{Act.F}, \code{area.region} and \code{area.module}.
#' The parameter \code{fitting} refers to Moore-Penrose (MP) or Moore-Penrose L2 penalty (MP+L2). "MP" is the default option.
#' The parameter \code{neurons} refers to the number of hidden neurons, default is 100. The paramter \code{Act.F} refers to the activation function of the hidden neurons being a sigmoidal neuron the default and only option: Act.F = 'sig'.
#' The parameters \code{area.region} and \code{area.module} are necessary if you want to downscale with a variant of ELM, called Receptive Fields ELM (RF-ELM). The parameter \code{area.region} is a vector with two parameters (i.e. c(a,b)), meaning the number of consecutive points of x in latitude (a) and in longitude (b). Default is NULL. The parameter \code{area.module}, is the size of the area within the grid region that is masked and fed to the hidden neurons. Default is NULL.
#' 
#' \strong{Neural Networks}
#' 
#' Neural network is based on the library \pkg{deepnet}. The optional parameters corresponds to those in \code{\link[deepnet]{nn.train}} and are: \code{initW} = NULL, \code{initB} = NULL, \code{hidden} = c(10), \code{activationfun} = "sigm", \code{learningrate} = 0.8, \code{momentum} = 0.5, \code{learningrate_scale} = 1, \code{output} = "sigm", \code{numepochs} = 3, \code{batchsize} = 100, \code{hidden_dropout} = 0, \code{visible_dropout} = 0. The values indicated are the default values.
#' 
#' \strong{Help}
#' 
#' If there are still doubts about the optional parameters despite the description here, we encourage to look for further details in the atomic functions: \code{\link[downscaleR]{glm.train}}, \code{\link[downscaleR]{mlr.train}}, \code{\link[downscaleR]{elm.train}} and \code{\link[deepnet]{nn.train}}.
#' 
#' @return A list of objects that contains the prediction on the train dataset, the model, the predictors and predictands used.
#' \itemize{
#'    \item \code{pred}: An object with the same structure as the predictands input parameter, but with pred$Data being the predictions and not the observations.
#'    \item \code{model}: A list with the information of the model: method, coefficients, fitting technique...
#'    \item \code{predictors}: Same as the predictors input parameter.
#'    \item \code{predictands}: Same as the predictands input parameter.}
#'    
#' @author J. Bano-Medina
#' @export
#' @importFrom MASS ginv
#' @importFrom matlab reshape repmat
#' @import deepnet 
#' @examples 
#' # Loading predictors
#' x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_psl, NCEP_Iberia_ta850)
#' x <- subsetGrid(x, years = 1985:1995)
#' # Loading predictands
#' y <- VALUE_Iberia_pr
#' # Prepare predictors
#' xtrain <- prepare_predictors(x = x, y = y)
#' # Downscaling a generalized linear model. In particular a logistic regression.
#' # model <- downscale.train(ybin, xtrain, method = "GLM" ,family = binomial(link = "logit"))
#' # Downscaling a linear model..
#' model <- downscale.train(xtrain$y, xtrain, method = "MLR", fitting = 'MP')
#' str(model)
#' # Downscaling a non-linear model via a regularized L2 Extreme learning Machine..
#' model <- downscale.train(xtrain$y, xtrain, method = "ELM", fitting = 'MP+L2')
#' str(model)
#' # Downscaling a non-linear model via Neural Networks..
#' model <- downscale.train(xtrain$y, xtrain, method = "NN", learningrate = 0.1, numepochs = 10, hidden = 5, output = 'linear') 
#' str(model)
#' # Downscaling a local model with the closest 4 grid points.
#' xtrain <- prepare_predictors(x = x,y = y,local.predictors = list(neigh.vars = "shum850",n.neighs = 4))
#' xtest <- prepare_newdata(newdata = x, predictor = xtrain)
#' model <- downscale.train(xtrain$y, xtrain, method = "MLR", fitting = 'MP')
#' str(model)
#' # Donwscaling with Principal Components (PCs)
#' xtrain <- prepare_predictors(x = x,y = y, PCA = list(v.exp = 0.9))
#' xtest <- prepare_newdata(newdata = x, predictor = xtrain)
#' model <- downscale.train(xtrain$y, xtrain, method = "ELM" ,fitting = 'MP')
#' str(model)
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

##############################################################################################################
#                     DOWNSCALING                                                                            #
##############################################################################################################
#' @title Switch to selected downscale method.
#' @description Internal function of \code{\link[downscaleR]{downscale.train}} that switches to the selected method.
#' @param x The input grid. Class: matrix.
#' @param y The observations dataset. Class: matrix.
#' @param method Type of transer function. Options are: GLM, MLR, ELM and NN. The default is MLR.
#' @param ... Optional parameters. These parameters are different depending on the method selected. Every parameter has a default value set in the atomic functions in case that no selection is wanted. For this reason see the atomic functions for more details: \code{\link[downscaleR]{glm.train}}, \code{\link[downscaleR]{mlr.train}}, \code{\link[downscaleR]{elm.train}} and \code{\link[deepnet]{nn.train}}.  
#' @return An object with the information of the selected model.
#' @details The optional parameters of neural networks can be found in the library \pkg{deepnet} via \code{\link[deepnet]{nn.train}}This function is internal and should not be used by the user. The user should use \code{\link[downscaleR]{downscale.train}}.
#' @author J. Bano-Medina
#' @export
downs.train <- function(x, y, method, ...) {
  switch(method,
         GLM   = atomic_model <- glm.train(x, y, ...),
         MLR   = atomic_model <- mlr.train(x, y, ...),
         ELM   = atomic_model <- elm.train(x, y, ...),
         NN    = atomic_model <- nn.train(x, y,  ...))
  return(atomic_model)}

##############################################################################################################
#                     Generalized Linear Models (GLM)                                                        #
##############################################################################################################
#' @title Donwscaling with generalized linear models (GLM).
#' @description Donwscaling with generalized linear models (GLM) with the base function \code{\link[stats]{glm}}.
#' @param x The grid data. Class: matrix.
#' @param y The observations data. Class: matrix.
#' @param family  Family of the generalized linear model according to \code{\link[stats]{family}}. This is an optional parameter and there is no default. 
#' @return The GLM model as returned from \code{\link[stats]{glm}}.
#' @details The indicated optional parameters MUST be included in \code{\link[downscaleR]{downscale.predict}}, in this case \code{family} is the only optional parameter. This function is internal and should not be used by the user. The user should use \code{\link[downscaleR]{downscale.predict}}.
#' @author J. Bano-Medina
#' @export
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

##############################################################################################################
#                     Multiple Linear Regression (GLM)                                                       #
##############################################################################################################
#' @title Donwscaling with multiple linear regression (MLR).
#' @description Donwscaling with multiple linear regression (MLR).
#' @param x The grid data. Class: matrix.
#' @param y The observations data. Class: matrix.
#' @param fitting  Type of fitting used: c("LS","MP","MP+L2"). There are three options: least squares ("LS") that downscales via \code{\link[stats]{lm}}, Moore-Penrose ("MP") that uses the generalized inverse via \code{\link[MASS]{ginv}} or the L2 penalty, commonly known as ridge regression ("MP+L2"). This is an optional parameter and the default is "MP". 
#' @return If fitting = "LS", then returns a MLR model as returned from \code{\link[stats]{lm}} whereas if fitting = ("MP" or "MP+L2") then returns just the estimated coefficients/weights.
#' @details If "MP" or "MP+L2" are selected, the bias is input as an additional input neuron in the weight matrix with value equal 1. The indicated optional parameters can be changed by including them in \code{\link[downscaleR]{downscale.predict}}, in this case \code{fitting} is the only optional parameter. This function is internal and should not be used by the user. The user should use \code{\link[downscaleR]{downscale.predict}}.
#' @author J. Bano-Medina
#' @export
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

##############################################################################################################
#                     Extreme Learning Machine (ELM)                                                         #
##############################################################################################################
#' @title Donwscaling with extreme learning machine (ELM).
#' @description Donwscaling with extreme learning machine (ELM).
#' @param x The grid data. Class: matrix.
#' @param y The observations data. Class: matrix.
#' @param fitting  Type of fitting used: c(MP","MP+L2"). There are two options: Moore-Penrose ("MP") that uses the generalized inverse via \code{\link[MASS]{ginv}} or the L2 penalty, commonly known as ridge regression ("MP+L2"). This is an optional parameter and the default is "MP". 
#' @param neurons The number of hidden neurons. This is an optional parameter and the default is 100. 
#' @param Act.F The type of activation function of the hidden neurons. This is an optional parametern and the default is 'sig' or sigmoidal.
#' @param area.region A vector with two parameters (i.e. c(a,b)), meaning the number of consecutive points of x in latitude (a) and in longitude (b). This is an optional parameter and the default is NULL.
#' @param area.module The size of the area within the grid region that is masked and fed to the hidden neurons. This is an optional parameter and the default is NULL.
#' @return If fitting = "LS", then returns a MLR model as returned from \code{\link[stats]{lm}} whereas if fitting = ("MP" or "MP+L2") then returns just the estimated coefficients/weights.
#' @details  The bias is input as an additional input neuron in the weight matrix with value equal 1. The indicated optional parameters can be changed by including them in \code{\link[downscaleR]{downscale.predict}}. This function is internal and should not be used by the user. The user should use \code{\link[downscaleR]{downscale.predict}}.
#' @author J. Bano-Medina

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