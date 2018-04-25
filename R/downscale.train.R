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
#' @description Downscale data to local scales by statistical methods: analogs, generalized linear models (GLM) and Neural Networks (NN). 
#' @param obj An object. The object as returned by \code{\link[downscaleR]{prepareData}}.
#' @param method A string value. Type of transer function. Options are c("analogs","GLM","NN").
#' @param filter A logical expression (i.e. = ">0"). This will filter all values that do not accomplish that logical statement. Default is NULL.
#' @param ... Optional parameters. These parameters are different depending on the method selected. Every parameter has a default value set in the atomic functions in case that no selection is wanted. 
#' Everything concerning these parameters is explained in the section \code{Details}. 
#' However, if wanted, the atomic functions can be seen here: \code{\link[downscaleR]{glm.train}} and \code{\link[deepnet]{nn.train}}.  

#' @details The function can downscale in both global and local mode, though not simultaneously. 
#' If there is perfect collinearity among predictors, then the matrix will not be invertible and the downscaling will fail.
#' We recommend to get rid of the NaN/NA values before calling the function.
#' 
#' \strong{Analogs}
#' The optional parameters of this method are:
#' \itemize{
#' \item \code{n.analogs} An integer. Number of analogs. Default is 4.
#' \item \code{sel.fun} A string. Select a function to apply to the analogs selected for a given observation. Options are 
#' "mean", "wmean" (i.e., weighted mean), "max", "min", "median", "prcXX" 
#' (i.e., prc85 means the 85th percentile of the analogs values distribution). Default is "mean".
#' the function applied to the analogs values, (i.e., sel.fun = c("mean","max","min","median","prcXX"), with default "mean") 
#' and the temporal window, (i.e., window = 0).
#' \item \code{window} An integer. Window of days removed when selecting analogs. 
#' If window = 7, then 7 days after the observation date and the 7 days before the observation date are removed. Default is 0.
#' \item \code{n.random} An integer. Choose N random analogs among the closest n.analogs. Default is NULL.
#' }
#' More information can be found in \code{\link[downscaleR]{analogs.train}}
#' 
#' \strong{Generalized Linear Models (GLM)}
#' The optional parameters depends on the \code{fitting} optional parameter:
#' \itemize{
  #' \item \code{fitting} A string indicating the types of objective functions and how to fit the linear model.
  #' \itemize{
    #' \item \code{fitting = NULL} In this case the generalized linear model uses the \code{\link[stats]{glm}} function to fit the linear model. 
    #' This is the default option.
    #' The optional parameters when fitting = NULL are:
    #' \itemize{
      #' \item \code{family} A string indicating a description of the error distribution. Options are 
      #' family = c("gaussian","binomial","Gamma","inverse.gaussian","poisson","quasi","quasibinomial","quasipoisson"). 
      #' The links can be also specified and can be found in \code{\link[stats]{family}}.
      #' \item \code{na.action} A function which indicates what should happen when the data contain NAs. 
      #' The default is set by the na.action setting of options, and is na.fail if that is unset. 
      #' The ‘factory-fresh’ default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
      #' }
    #' \item \code{fitting = "stepwise"} Indicates a stepwise regression via \code{\link[stats]{glm}} and \code{\link[stats]{step}}.
    #' The optional parameters are the same than for fitting = NULL. The stepwise performs always a forward selection search stopping
    #' \item \code{fitting = c("L1","L2","L1L2","gLASSO")}. These four options refer to ridge regression (L1 penalty), lasso regression (L2 penalty),
    #' elastic-net regression (L1L2 penalty) and group Lasso regression (group L2 penalty). The model is fitted via 
    #' \code{\link[glmnet]{glmnet}} and the corresponding penalties are found via \code{\link[glmnet]{cv.glmnet}}. This function \code{\link[glmnet]{glmnet}}
    #' forces by default to standardize predictors, however we have changed it to standardize = FALSE, and standardization should be done prior to 
    #' the downscaling process. 
    #' The optional parameters when fitting = c("L1","L2","L1L2","gLASSO") are:
    #' \itemize{
      #' \item \code{family} A string indicating a description of the error distribution. Options are 
      #' family = c("gaussian","binomial","Gamma","inverse.gaussian","poisson","quasi","quasibinomial","quasipoisson"). 
      #' The links CAN NOT be specified as the \code{\link[glmnet]{glmnet}} has not been programmed to handle links.
      #' However, the default ones can be found in \code{\link[stats]{family}}. If fitting = "gLASSO" then family must be "mgaussian".
      #' \item \code{offset} A vector of length nobs that is included in the linear predictor (a nobs x nc matrix for the "multinomial" family). 
      #' Useful for the "poisson" family (e.g. log of exposure time), or for refining a model by starting at a current fit. 
      #' Default is NULL. If supplied, then values must also be supplied to the predict function.
      #' }
   #' }
#' There are two things to consider. 
#' 1) If family = "binomial" then type = "response" when predicting values.
#' 2) Except for fitting = "MP", for the rest of the fitting options, the parameter site must be TRUE, unless 
#' we want a gLASSO, in this case site must be FALSE.
#' 
#' }
#' 
#' \strong{Neural Networks}
#' Neural network is based on the library \pkg{deepnet}. The optional parameters corresponds to those in \code{\link[deepnet]{nn.train}}
#' and are: \code{initW} = NULL, \code{initB} = NULL, \code{hidden} = c(10), \code{activationfun} = "sigm", \code{learningrate} = 0.8, \code{momentum} = 0.5, 
#' \code{learningrate_scale} = 1, \code{output} = "sigm", \code{numepochs} = 3, \code{batchsize} = 100, \code{hidden_dropout} = 0, \code{visible_dropout} = 0. The values indicated are the default values.
#' 
#' \strong{Help}
#' 
#' If there are still doubts about the optional parameters despite the description here, we encourage to look for further details in the atomic functions: 
#' \code{\link[downscaleR]{analogs.train}}, \code{\link[downscaleR]{glm.train}} and \code{\link[deepnet]{nn.train}}.
#' 
#' @return A list of objects that contains the prediction on the train dataset and the model.
#' \itemize{
#'    \item \code{pred}: An object with the same structure as the predictands input parameter, but with pred$Data being the predictions and not the observations.
#'    \item \code{model}: A list with the information of the model: method, coefficients, fitting ...
#'    }
#' @importFrom transformeR isRegular
#' @author J. Bano-Medina
#' @export
#' @importFrom MASS ginv
#' @import deepnet 
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
#' # Downscaling PRECIPITATION
#' # ... via analogs ...
#' model <- downscale.train(xyT, method = "analogs",
#'                          sel.fun = "mean")
#' # ... via a logistic regression (ocurrence of precipitation)
#' # and gaussian regression (amount of precipitation) ...
#' model.ocu <- downscale.train(xyT.bin, method = "GLM",
#'                              family = binomial(link = "logit"))
#' model.reg <- downscale.train(xyT, method = "GLM",
#'                              family = "gaussian", filter = ">0")
#' # ... via a neural network ...
#' model.ocu <- downscale.train(xyT.bin, method = "NN",
#'                              learningrate = 0.1, numepochs = 10, hidden = 5,
#'                              output = 'linear')
#' model.reg <- downscale.train(xyT, method = "NN",
#'                              learningrate = 0.1, numepochs = 10,
#'                              hidden = 5, output = 'linear')
#' # Downscaling PRECIPITATION - Local model with the closest
#' # 4 grid points and multisite linear regression.
#' xyT.local <- prepareData(x = x, y = y,
#'                          local.predictors = list(vars = "hus@850",n = 4))
#' model <- downscale.train(xyT.local,method = "analogs")
#' # Downscaling PRECIPITATION - Principal Components (PCs)
#' # and gamma regression for the amount of precipitation
#' xyT.pc     <- prepareData(x = x,y = y,
#'                           spatial.predictors = list(which.combine = getVarNames(x),v.exp = 0.9))
#' model <- downscale.train(xyT.local,method = "analogs")

downscale.train <- function(obj, method, filter = NULL, ...) {
  dimNames <- getDim(obj$y)
  pred <- obj$y
  if ( method == "GLM") {
    if (attr(obj,"nature") == "spatial+local") {
      site <- "mix"}
    else {
      site <- "single"
    }
  }
  if ( method == "NN" || method == "analogs") {
    if (attr(obj,"nature") == "local") {
      site <- "single"
    }
    else if (attr(obj,"nature") == "spatial+local") {
      site <- "mix"
    }
    else {
      site <- "multi"
    }
  }
  
  # Multi-site
  if (site == "multi") {
    yy <- obj$y$Data
    if (method == "analogs") {
      atomic_model <- downs.train(obj$x.global, yy, method, dates = getRefDates(obj$y), ...)}
    else {      
      atomic_model <- downs.train(obj$x.global, yy, method, ...)}
    if (method == "analogs") {atomic_model$dates$test <- getRefDates(obj$y)}
    pred$Data <- downs.predict(obj$x.global, method, atomic_model)}
  # Single-site
  else if (site == "single") {
    pred$Data    <- array(data = NA, dim = dim(obj$y$Data)); attr(pred$Data,"dimensions") <- attr(obj$y$Data,"dimensions")
    regular <- FALSE
    if (isRegular(obj$y)) {
      pred$Data <- array3Dto2Dmat(pred$Data)
      obj$y$Data <- array3Dto2Dmat(obj$y$Data)
      regular <- TRUE
    }
    stations <- dim(pred$Data)[2]
    atomic_model <- vector("list",stations)
    for (i in 1:stations) {
      if (attr(obj,"nature") == "local") {
        xx = obj$x.local[[i]]$member_1}
      else {
        xx = obj$x.global}
      yy = obj$y$Data[,i, drop = FALSE]
      if (is.null(filter)) {ind = eval(parse(text = "which(!is.na(yy))"))}
      else {ind = eval(parse(text = paste0("which(!is.na(yy) & yy",filter,")")))}
      if (method == "analogs") {
        atomic_model[[i]] <- downs.train(xx[ind,, drop = FALSE], yy[ind,,drop = FALSE], method, dates = getRefDates(obj$y)[ind], ...)}
      else {
        atomic_model[[i]] <- downs.train(xx[ind,, drop = FALSE], yy[ind,,drop = FALSE], method, ...)}
      if (method == "analogs") {atomic_model[[i]]$dates$test <- getRefDates(obj$y)}
      pred$Data[,i] <- downs.predict(xx, method, atomic_model[[i]])
    }
    if (regular) {
      pred$Data <- mat2Dto3Darray(pred$Data, x = pred$xyCoords$x, y = pred$xyCoords$y)
    }
  }
  # Mix - Global predictors with local predictors
  else if (site == "mix") {
    pred$Data    <- array(data = NA, dim = dim(obj$y$Data)); attr(pred$Data,"dimensions") <- attr(obj$y$Data,"dimensions")
    regular <- FALSE
    if (isRegular(obj$y)) {
      pred$Data <- array3Dto2Dmat(pred$Data)
      obj$y$Data <- array3Dto2Dmat(obj$y$Data)
      regular <- TRUE
    }
    stations <- dim(pred$Data)[which(getDim(pred) == "loc")]
    atomic_model <- vector("list",stations)
    for (i in 1:stations) {
      xx1 = obj$x.local[[i]]$member_1
      xx2 = obj$x.global
      xx <- cbind(xx1,xx2)
      yy = obj$y$Data[,i, drop = FALSE]
      if (is.null(filter)) {ind = eval(parse(text = "which(!is.na(yy))"))}
      else {ind = eval(parse(text = paste0("which(!is.na(yy) & yy",filter,")")))}
      if (method == "analogs") {
        atomic_model[[i]] <- downs.train(xx[ind,, drop = FALSE], yy[ind,,drop = FALSE], method, dates = getRefDates(obj$y)[ind], ...)}
      else {
        atomic_model[[i]] <- downs.train(xx[ind,, drop = FALSE], yy[ind,,drop = FALSE], method, ...)}
      if (method == "analogs") {atomic_model[[i]]$dates$test <- getRefDates(obj$y)}
      pred$Data[,i] <- downs.predict(xx, method, atomic_model[[i]])}
    if (regular) {
      pred$Data <- mat2Dto3Darray(pred$Data, x = pred$xyCoords$x, y = pred$xyCoords$y)
    }
  }
  
  attr(pred$Data, "dimensions") <- dimNames
  model <- list("pred" = pred, "model" = list("method" = method, "site" = site, "atomic_model" = atomic_model))
  return(model)}

##############################################################################################################
#                     DOWNSCALING                                                                            #
##############################################################################################################
#' @title Switch to selected downscale method.
#' @description Internal function of \code{\link[downscaleR]{downscale.train}} that switches to the selected method.
#' @param x The input grid. Class: matrix.
#' @param y The observations dataset. Class: matrix.
#' @param method Type of transer function. Options are: analogs, GLM and NN. 
#' @param ... Optional parameters. These parameters are different depending on the method selected. Every parameter has a default value set in the atomic functions in case that no selection is wanted. For this reason see the atomic functions for more details: \code{\link[downscaleR]{glm.train}} and \code{\link[deepnet]{nn.train}}.  
#' @return An object with the information of the selected model.
#' @details The optional parameters of neural networks can be found in the library \pkg{deepnet} via \code{\link[deepnet]{nn.train}}This function is internal and should not be used by the user. The user should use \code{\link[downscaleR]{downscale.train}}.
#' @author J. Bano-Medina
#' @export
downs.train <- function(x, y, method, ...) {
  switch(method,
         analogs = atomic_model <- analogs.train(x, y, ...),
         GLM     = atomic_model <- glm.train(x, y, ...),
         NN      = atomic_model <- nn.train(x, y,  ...))
  return(atomic_model)}

