##     downscaleCV.R Downscaling method calibration in cross validation mode.
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

#' @title Downscale climate data and reconstruct the temporal serie by splitting the data following a user-defined scheme
#' @description Downscale climate data and reconstruct the temporal serie by splitting the data following a user-defined scheme.
#' The statistical downscaling methods currently implemented are: analogs, generalized linear models (GLM) and Neural Networks (NN). 
#' @param x The input grid (admits both single and multigrid, see \code{\link[transformeR]{makeMultiGrid}}). It should be an object as returned by \pkg{loadeR}.
#' @param y The observations dataset. It should be an object as returned by \pkg{loadeR}.
#' @param method A string value. Type of transer function. Currently implemented options are \code{"analogs"}, \code{"GLM"} and \code{"NN"}.
#' @param simulate A logic value indicating whether we want to simulate or not based on the GLM distributional parameters. Only relevant when perdicting with a GLM. Default to FALSE. 
#' @param sampling.strategy Specifies a sampling strategy to define the training and test subsets. Possible values are 
#' \code{"kfold.chronological"} (the default), \code{"kfold.random"}, \code{"leave-one-year-out"} and NULL.
#' The \code{sampling.strategy} choices are next described:
#' \itemize{
#'   \item \code{"kfold.random"} creates the number of folds indicated in the \code{folds} argument by randomly sampling the entries along the time dimension.
#'   \item \code{"kfold.chronological"} is similar to \code{"kfold.random"}, but the sampling is performed in ascending order along the time dimension.
#'   \item \code{"leave-one-year-out"}. This scheme performs a leave-one-year-out cross validation. It is equivalent to introduce in the argument \code{folds} a list of all years one by one.
#'   \item \code{NULL}. The folds are specified by the user in the function parameter \code{folds}.
#' }
#' The first two choices will be controlled by the argument \code{folds} (see below)
#' @param folds This arguments controls the number of folds, or how these folds are created (ignored if \code{sampling.strategy = "leave-one-year-out"}). If it is given as a fraction in the range (0-1), 
#' it splits the data in two subsets, one for training and one for testing, being the given value the fraction of the data used for training (i.e., 0.75 will split the data so that 75\% of the instances are used for training, and the remaining 25\% for testing). 
#' In case it is an integer value (the default, which is 4), it sets the number of folds in which the data will be split (e.g., \code{folds = 10} for the classical 10-fold cross validation). 
#' Alternatively, this argument can be passed as a list, each element of the list being a vector of years to be included in each fold (See examples).
#' @param scaleGrid.args A list of the parameters related to scale grids. This parameter calls the function \code{\link[transformeR]{scaleGrid}}. See the function definition for details on the parameters accepted.
#' @param prepareData.args A list with the arguments of the \code{\link[downscaleR]{prepareData}} function. Please refer to \code{\link[downscaleR]{prepareData}} help for
#' more details about this parameter.
#' @param condition Inequality operator to be applied considering the given threshold.
#' \code{"GT"} = greater than the value of \code{threshold}, \code{"GE"} = greater or equal,
#' \code{"LT"} = lower than, \code{"LE"} = lower or equal than. We only train with the days that satisfy the condition.
#' @param threshold Numeric value. Threshold used as reference for the condition. Default is NULL. If a threshold value is supplied with no specification of the parameter \code{condition}. Then condition is set to \code{"GE"}.
#' @param ... Optional parameters. These parameters are different depending on the method selected. 
#' Every parameter has a default value set in the atomic functions in case that no selection is wanted. 
#' Everything concerning these parameters is explained in the section \code{Details} of the function \code{\link[downscaleR]{downscaleTrain}}. However, if wanted, the atomic functions can be seen here: 
#' \code{\link[downscaleR]{analogs.train}}, \code{\link[downscaleR]{glm.train}} and \code{\link[deepnet]{nn.train}}.  
#' @details The function relies on \code{\link[downscaleR]{prepareData}}, \code{\link[downscaleR]{prepareNewData}}, \code{\link[downscaleR]{downscaleTrain}}, and \code{\link[downscaleR]{downscalePredict}}. 
#' For more information please visit these functions. It is envisaged to allow for a flexible fine-tuning of the cross-validation scheme. It uses internally the \pkg{transformeR} 
#' helper \code{\link[transformeR]{dataSplit}} for flexible data folding. 
#' Note that the indices for data splitting are obtained using \code{\link[transformeR]{getYearsAsINDEX}} when needed (e.g. in leave-one-year-out cross validation), 
#' thus adequately handling potential inconsistencies in year selection when dealing with year-crossing seasons (e.g. DJF).
#' 
#' If the variable to downscale is the precipitation and it is a binary variable,
#'  then two temporal series will be returned:
#' \enumerate{
#' \item The temporal serie with binary values filtered by a threshold adjusted by the train dataset, see \code{\link[transformeR]{binaryGrid}} for more details.
#' \item The temporal serie with the results obtained by the downscaling, without binary conversion process.
#' }
#' 
#' Missing data removal is recommended prior to multisite calibration.
#' 
#' According to the concept of cross-validation, a particular year should not appear in more than one fold
#' (i.e., folds should constitute disjoint sets). For example, the choice \code{fold =list(c(1988,1989), c(1989, 1990))}
#'  will raise an error, as 1989 appears in more than one fold.
#' 
#' @return The reconstructed downscaled temporal serie.
#' @seealso 
#' downscaleTrain for training a downscaling model
#' downscalePredict for prediction for a a test dataset with a trained model for 
#' \href{https://github.com/SantanderMetGroup/downscaleR/wiki/training-downscaling-models}{downscaleR Wiki} for downscaling seasonal forecasting and climate projections.
#' @importFrom transformeR dataSplit scaleGrid binaryGrid makeMultiGrid
#' @family downscaling.functions
#' @author J. Bano-Medina
#' @export
#' @examples \donttest{
#' require(transformeR)
#' require(climate4R.datasets)
#' data(NCEP_Iberia_hus850, NCEP_Iberia_ta850)
#' x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_ta850)
#' x <- subsetGrid(x, years = 1985:1995)
#' # Loading predictands
#' data("VALUE_Iberia_pr")
#' y <- VALUE_Iberia_pr
#' y <- getTemporalIntersection(obs = y, prd = x, "obs")
#' x <- getTemporalIntersection(obs = y, prd = x, "prd")
#' # ... kfold in 3 parts equally divided ...
#' pred <- downscaleCV(x, y, folds = 3, sampling.strategy = "kfold.chronological",
#'                      scaleGrid.args = list(type = "standardize"),
#'                      method = "GLM", family = Gamma(link = "log"), condition = "GT", threshold = 0,
#'                      prepareData.args = list(
#'                      "spatial.predictors" = list(which.combine = getVarNames(x), v.exp = 0.9)))
#' # ... kfold by years ...
#' pred <- downscaleCV(x,y,sampling.strategy = "kfold.chronological",
#'                      method = "GLM", condition = "GT", threshold = 0,
#'                      scaleGrid.args = list(type = "standardize"),
#'                      folds = list(c(1985,1986,1987,1988),
#'                                   c(1989,1990,1991,1992),
#'                                   c(1993,1994,1995)))
#' # ... leave one year out  ...
#' pred <- downscaleCV(x,y,sampling.strategy = "leave-one-year-out",
#'                      method = "GLM", condition = "GT", threshold = 0,
#'                      scaleGrid.args = list(type = "standardize"))
#' # Reconstructing the downscaled serie in 3 folds with local predictors.
#' pred <- downscaleCV(x,y,folds = 3, sampling.strategy = "kfold.chronological",
#'                      scaleGrid.args = list(type = "standardize"),
#'                      method = "GLM", condition = "GT", threshold = 0,
#'                      prepareData.args = list("local.predictors" = list(vars = "hus@850", n = 4)))
#' }

downscaleCV <- function(x, y, method,
                        sampling.strategy = "kfold.chronological", folds = 4, 
                        scaleGrid.args = NULL, simulate = FALSE,
                        prepareData.args = list("global.vars" = NULL, "combined.only" = TRUE, "spatial.predictors" = NULL, "local.predictors" = NULL, "extended.predictors" = NULL),
                        condition = NULL, threshold = NULL, ...) {
  
  if (!exists("global.vars",prepareData.args)) prepareData.args$global.vars <- NULL
  if (!exists("combined.only",prepareData.args)) prepareData.args$combined.only <- TRUE
  if (!exists("spatial.predictors",prepareData.args)) prepareData.args$spatial.predictors <- NULL
  if (!exists("local.predictors",prepareData.args)) prepareData.args$local.predictors <- NULL
  if (!exists("extended.predictors",prepareData.args)) prepareData.args$extended.predictors <- NULL
  
  x <- getTemporalIntersection(x,y,which.return = "obs")
  y <- getTemporalIntersection(x,y,which.return = "prd")
  
  if (!is.null(sampling.strategy)) {
    if (sampling.strategy == "leave-one-year-out") {
      type <- "chronological"
      folds <- as.list(getYearsAsINDEX(y) %>% unique())
    }
    
    if (sampling.strategy == "kfold.chronological") {
      type <- "chronological"
      if (!is.numeric(folds)) {
        folds.user <- unlist(folds) %>% unique() %>% sort()
        folds.data <- getYearsAsINDEX(y) %>% unique()
        if (any(folds.user != folds.data)) stop("In the parameters folds you have indicated years that do not belong to the dataset. Please revise the setup of this parameter.")
      }
    }
    if (sampling.strategy == "kfold.random") {
      type <- "random"
      if (!is.numeric(folds)) stop("In kfold.random, the parameter folds represent the NUMBER of folds and thus, it should be a NUMERIC value.")
    }
  }
  if (is.list(folds)) {
    if (any(duplicated(unlist(folds)))) stop("Years can not appear in more than one fold")
  }
  
  data <- dataSplit(x,y, f = folds, type = type)
  p <- lapply(1:length(data), FUN = function(xx) {
    message(paste("fold:",xx,"-->","calculating..."))
    xT <- data[[xx]]$train$x ; yT <- data[[xx]]$train$y
    xt <- data[[xx]]$test$x  ; yt <- data[[xx]]$test$y
    if (!is.null(scaleGrid.args)) {
      scaleGrid.args$base <- xT
      scaleGrid.args$grid <- xt
      scaleGrid.args$skip.season.check <- TRUE
      xt <- do.call("scaleGrid",args = scaleGrid.args)
      scaleGrid.args$grid <- xT
      xT <- do.call("scaleGrid",args = scaleGrid.args)
    }
    xT <- prepareData(x = xT, y = yT,  global.vars = prepareData.args$global.vars, combined.only = prepareData.args$combined.only, spatial.predictors = prepareData.args$spatial.predictors,local.predictors = prepareData.args$local.predictors,extended.predictors = prepareData.args$extended.predictors)
    xt <- prepareNewData(newdata = xt, data.structure = xT)
    model <- downscaleTrain(xT, method, condition, threshold, ...)
    if (all(as.vector(y$Data) %in% c(0,1,NA,NaN), na.rm = TRUE)) {
      y.prob <- downscalePredict(xt, model, simulate)
      if (method == "GLM") {
        if (isTRUE(simulate)) {
          y.bin  <- y.prob
        }
        else {
          y.bin  <- binaryGrid(y.prob, ref.obs = yT, ref.pred = model$pred)
        }
      }
      else{
        y.bin  <- binaryGrid(y.prob, ref.obs = yT, ref.pred = model$pred)}
      out <- makeMultiGrid(list(y.prob,y.bin)) %>% redim(drop = TRUE)
      out$Variable$varName <- c("prob","bin")
    }
    else {
      out <- downscalePredict(xt, model, simulate)
    }
    return(out)
  }) %>% bindGrid(dimension = "time")
  return(p)
}
