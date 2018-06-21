##     downscale.cv.R Downscaling method calibration in cross validation mode.
##
##     Copyright (C) 2018 Santander Meteorology Group (http://www.meteo.unican.es)
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

#' @title Downscale climate data and reconstruct the temporal serie by splitting the data in k folds.
#' @description Downscale climate data and reconstruct the temporal serie by splitting the data in k folds. 
#' Statistical downscaling methods are: analogs, generalized linear models (GLM) and Neural Networks (NN). 
#' @param x The input grid (admits both single and multigrid, see \code{\link[transformeR]{makeMultiGrid}}). It should be an object as returned by \pkg{loadeR}.
#' @param y The observations dataset. It should be an object as returned by \pkg{loadeR}.
#' @param method A string value. Type of transer function. Currently implemented options are \code{"analogs"}, \code{"GLM"} and \code{"NN"}.
#' @param sampling.strategy Specifies a sampling strategy to define the training and test subsets. Possible values are \code{"kfold.random"}, \code{"kfold.chronological"} and \code{"leave-one-year-out"}.
#' The \code{sampling.strategy} choices are next described:
#' \itemize{
#'   \item \code{"kfold.random"} creates the number of folds indicated in the \code{folds} argument by randomly sampling the entries along the time dimension.
#'   \item \code{"kfold.chronological"} is similar to \code{"kfold.random"}, but the sampling is performed in ascending order along the time dimension.
#'   \item \code{"leave-one-year-out"}. This schema performs a leave-one-year-out cross validation. It is equivalent to introduce in the argument \code{folds} a list of all years one by one.
#' }
#' The first two choices will be controlled by the argument \code{folds} (see below)
#' @param folds This arguments controls the number of folds, or how these folds are created (ignored if \code{sampling.strategy = "leave-one-year-out"}). If it is given as a fraction in the range (0-1), 
#' it splits the data in two subsets, one for training and one for testing, being the given value the fraction of the data used for training (i.e., 0.75 will split the data so that 75\% of the instances are used for training, and the remaining 25\% for testing). 
#' In case it is an integer value, it sets the number of folds in which the data will be split (e.g., \code{folds = 10} for the classical 10-fold cross validation). 
#' Alternatively, this argument can be passed as a list, each element of the list being a vector of years to be included in each fold (See examples).
#' @param scale.list A list of the parameters related to scale grids. This parameter calls the function \code{\link[transformeR]{scaleGrid}}. See the function definition for details on the parameters accepted.
#' @param global.vars An optional character vector with the short names of the variables of the input x multigrid to be retained as global predictors 
#' (use the \code{\link[transformeR]{getVarNames}} helper if not sure about variable names). 
#' This argument just produces a call to \code{\link[transformeR]{subsetGrid}}, but it is included here for better flexibility in downscaling experiments (predictor screening...). 
#' For instance, it allows to use some specific variables contained in x as local predictors and the remaining ones, specified in subset.vars, 
#' as either raw global predictors or to construct the combined PC.
#' @param combined.only Optional, and only used if \code{spatial.predictors} parameters are passed. Should the combined PC be used as the only global predictor? Default to TRUE. 
#' Otherwise, the combined PC constructed with \code{which.combine} argument in \code{\link[transformeR]{prinComp}} is append to the PCs of the remaining variables within the grid.
#' @param spatial.predictors Default to \code{NULL}, and not used. Otherwise, a named list of arguments in the form argument = value, 
#' with the arguments to be passed to prinComp to perform Principal Component Analysis of the predictors grid (x). 
#' See Details on principal component analysis of predictors.
#' @param local.predictors Default to \code{NULL}, and not used. Otherwise, a named list of arguments in the form \code{argument = value},
#'  with the following arguments:
#'  \itemize{
#'    \item \code{vars}: names of the variables in \code{x} to be used as local predictors
#'    \item \code{fun}: Optional. Aggregation function for the selected local neighbours.
#'    The aggregation function is specified as a list, indicating the name of the aggregation function in
#'     first place (as character), and other optional arguments to be passed to the aggregation function.
#'     For instance, to compute the average skipping missing values: \code{fun = list(FUN= "mean", na.rm = TRUE)}.
#'     Default to NULL, meaning that no aggregation is performed.
#'    \item \code{n}: Integer. Number of nearest neighbours to use. If a single value is introduced, and there is more
#'    than one variable in \code{vars}, the same value is used for all variables. Otherwise, this should be a vector of the same
#'    length as \code{vars} to indicate a different number of nearest neighbours for different variables.
#'  }
#' @param extended.predictors This is a parameter related to the extreme learning machine and reservoir computing framework where input data is randomly projected into a new space of size \code{n}. Default to \code{NULL}, and not used. Otherwise, a named list of arguments in the form \code{argument = value},
#'  with the following arguments:
#'  \itemize{
#'    \item \code{n}: A numeric value. Indicates the size of the random nonlinear dimension where the input data is projected.
#'    \item \code{module}: A numeric value (Optional). Indicates the size of the mask's module. Belongs to a specific type of ELM called RF-ELM.
#'  }
#' @param condition Inequality operator to be applied considering the given threshold.
#' \code{"GT"} = greater than the value of \code{threshold}, \code{"GE"} = greater or equal,
#' \code{"LT"} = lower than, \code{"LE"} = lower or equal than. We only train with the days that satisfy the condition.
#' @param threshold Numeric value. Threshold used as reference for the condition. Default is NULL. If a threshold value is supplied with no specification of the parameter \code{condition}. Then condition is set to \code{"GE"}.
#' @param ... Optional parameters. These parameters are different depending on the method selected. 
#' Every parameter has a default value set in the atomic functions in case that no selection is wanted. 
#' Everything concerning these parameters is explained in the section \code{Details} of the function \code{\link[downscaleR]{downscale.train}}. However, if wanted, the atomic functions can be seen here: 
#' \code{\link[downscaleR]{glm.train}} and \code{\link[deepnet]{nn.train}}.  
#' @details The function relies on \code{\link[downscaleR]{prepareData}}, \code{\link[downscaleR]{prepareNewData}}, \code{\link[downscaleR]{downscale.train}}, and \code{\link[downscaleR]{downscale.predict}}. 
#' For more information please visit these functions.#' The function is envisaged to be flexible enough for allowing a fine-tuning of the cross-validation scheme. It uses internally the \pkg{transformeR} 
#' helper \code{\link[transformeR]{dataSplit}} for flexible data folding. 
#' Note that the indices for data splitting are obtained using \code{\link[transformeR]{getYearsAsINDEX}} when needed (e.g. in leave-one-year-out cross validation), 
#' thus adequately handling potential inconsistencies in year selection when dealing with year-crossing seasons (e.g. DJF).
#' If the variable to downscale is the precipitation and it is a binary variable, then two temporal series will be returned:
#' 1) The temporal serie with binary values filtered by a threshold adjusted by the train dataset, see \code{\link[transformeR]{binaryGrid}} for more details.
#' 2) The temporal serie with the results obtained by the downscaling, without any binary converting process.
#' We recommend to get remove missing data prior to multisite calibration.
#' According to the concept of cross-validation, a particular year should not appear in more than one fold. For example,
#' if fold.1 = c("1988","1989") and fold.2 = c("1989","1990") this will cause an error.
#' @return The reconstructed downscaled temporal serie.
#' @seealso 
#' downscale.train for training a downscaling model
#' downscale.predict for prediction for a a test dataset with a trained model for 
#' \href{https://github.com/SantanderMetGroup/downscaleR/wiki/training-downscaling-models}{downscaleR Wiki} for downscaling seasonal forecasting and climate projections.
#' @importFrom transformeR dataSplit scaleGrid binaryGrid
#' @family downscaling.functions
#' @author J. Bano-Medina
#' @export
#' @examples 
#' require(transformeR)
#' data(NCEP_Iberia_hus850, NCEP_Iberia_ta850)
#' x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_ta850)
#' x <- subsetGrid(x, years = 1985:1995)
#' # Loading predictands
#' y <- VALUE_Iberia_pr
#' y <- getTemporalIntersection(obs = y, prd = x, "obs")
#' x <- getTemporalIntersection(obs = y, prd = x, "prd")
#' # ... kfold in 3 parts equally divided ...
#' pred <- downscale.cv(x, y, folds = 3, sampling.strategy = "kfold.chronological",
#'                      scale.list = list(type = "standardize"),
#'                      method = "GLM", family = Gamma(link = "log"), condition = "GT", threshold = 0,
#'                      spatial.predictors = list(which.combine = getVarNames(x), v.exp = 0.9))
#' # ... kfold by years ...
#' pred <- downscale.cv(x,y,sampling.strategy = "kfold.chronological",
#'                      method = "GLM", condition = "GT", threshold = 0,
#'                      scale.list = list(type = "standardize"),
#'                      folds = list(c(1985,1986,1987,1988),
#'                                   c(1989,1990,1991,1992),
#'                                   c(1993,1994,1995)))
#' # ... leave one year out  ...
#' pred <- downscale.cv(x,y,sampling.strategy = "leave-one-year-out",
#'                      method = "GLM", condition = "GT", threshold = 0,
#'                      scale.list = list(type = "standardize"))
#' # Reconstructing the downscaled serie in 3 folds with local predictors.
#' pred <- downscale.cv(x,y,folds = 3, sampling.strategy = "kfold.chronological",
#'                      scale.list = list(type = "standardize"),
#'                      method = "GLM", condition = "GT", threshold = 0,
#'                      local.predictors = list(vars = "hus@850", n = 4))


downscale.cv <- function(x, y, method,
                         sampling.strategy = "kfold.chronological", folds = 4, 
                         scale.list = NULL,
                         global.vars = NULL, combined.only = TRUE, spatial.predictors = NULL, local.predictors = NULL, extended.predictors = NULL,
                         condition = NULL, threshold = NULL, ...) {
  
  x <- getTemporalIntersection(x,y,which.return = "obs")
  y <- getTemporalIntersection(x,y,which.return = "prd")
  
  if (sampling.strategy == "leave-one-year-out") {
    type <- "chronological"
    folds <- as.list(getYearsAsINDEX(y) %>% unique())
  }
  
  if (sampling.strategy == "kfold.chronological") {
    type <- "chronological"
    if (!is.numeric(folds)) {
      folds.user <- unlist(folds) %>% unique()
      folds.data <- getYearsAsINDEX(y) %>% unique()
      if (any(folds.user != folds.data)) stop("In the parameters folds you have indicated years that do not belong to the dataset. Please revise the setup of this parameter.")
    }
  }
  if (sampling.strategy == "kfold.random") {
    type <- "random"
    if (!is.numeric(folds)) stop("In kfold.random, the parameter folds represent the NUMBER of folds and thus, it should be a NUMERIC value.")
  }
  
  if (is.list(folds)) {
    if (any(duplicated(unlist(folds)))) stop("Years can not appear in more than one fold")
  }
  
  data <- dataSplit(x,y, f = folds, type = type)
  p <- lapply(1:length(data), FUN = function(xx) {
    message(paste("fold:",xx,"-->","calculating..."))
    xT <- data[[xx]]$train$x ; yT <- data[[xx]]$train$y
    xt <- data[[xx]]$test$x  ; yt <- data[[xx]]$test$y
    if (!is.null(scale.list)) {
      scale.list$base <- xT
      scale.list$grid <- xt
      scale.list$skip.season.check <- TRUE
      xt <- do.call("scaleGrid",args = scale.list)
      scale.list$grid <- xT
      xT <- do.call("scaleGrid",args = scale.list)
    }
    xT <- prepareData(x = xT, y = yT, global.vars = global.vars, combined.only = combined.only, spatial.predictors = spatial.predictors, local.predictors = local.predictors, extended.predictors = extended.predictors)
    xt <- prepareNewData(newdata = xt, data.structure = xT)
    model <- downscale.train(xT, method, condition, threshold, ...)
    if (all(as.vector(y$Data) %in% c(0,1,NA,NaN), na.rm = TRUE)) {
      y.prob <- downscale.predict(xt, model)
      if (method == "GLM") {
        if (model$model$atomic_model[[1]]$info$simulate == "yes") {
          y.bin  <- y.prob}
        else {
          y.bin  <- binaryGrid(y.prob, ref.obs = yT, ref.pred = model$pred)}}
      else{
        y.bin  <- binaryGrid(y.prob, ref.obs = yT, ref.pred = model$pred)}
      out <- list(y.prob$Data, y.bin$Data)}
    else{
      out <- list(downscale.predict(xt, model)$Data)}
    return(out)
    })
  
    pred <- lapply(1:length(p[[1]]), function(i) {
      pred <- y
      dimNames <- getDim(pred)
      pp <- lapply(1:length(data), function(z) p[[z]][[i]])  ; pred$Data   <- do.call(rbind,as.matrix(pp))
      attr(pred$Data, "dimensions") <- dimNames  
      return(pred)})
    if (length(pred) == 1) pred <- pred[[1]]
  return(pred)}
