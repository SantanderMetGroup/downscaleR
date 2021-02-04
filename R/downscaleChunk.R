##     downscaleChunk.R Downscaling method calibration in cross validation mode.
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

#' @title Downscaling by chunks
#' @description Downscale climate data by splitting it in chunks, where there are as many chunks as latitudes. 
#' This function encapsulates \code{\link[downscaleR]{downscaleTrain}}
#' @param x The input grid (admits both single and multigrid, see \code{\link[transformeR]{makeMultiGrid}}). It should be an object as returned by \pkg{loadeR}.
#' @param y The observations dataset. It should be an object as returned by \pkg{loadeR}.
#' @param newdata New datasets where to apply the model infered. It should be a list of objects as returned by \pkg{loadeR},
#' containing the new dataset/s.
#' @param method A string value. Type of transer function. Currently implemented options are \code{"analogs"}, \code{"GLM"} and \code{"NN"}.
#' @param simulate A logic value indicating whether we want to simulate or not based on the GLM distributional parameters. Only relevant when perdicting with a GLM. Default to FALSE. 
#' @param prepareData.args A list with the arguments of the \code{\link[downscaleR]{prepareData}} function. Please refer to \code{\link[downscaleR]{prepareData}} help for
#' more details about this parameter.
#' @param condition Inequality operator to be applied considering the given threshold.
#' \code{"GT"} = greater than the value of \code{threshold}, \code{"GE"} = greater or equal,
#' \code{"LT"} = lower than, \code{"LE"} = lower or equal than. We only train with the days that satisfy the condition.
#' @param threshold Numeric value. Threshold used as reference for the condition. Default is NULL. If a threshold value is supplied with no specification of the parameter \code{condition}. Then condition is set to \code{"GE"}.
#' @param predict A logic value. Should the prediction on the training set should be returned? Default is TRUE.
#' @param ... Optional parameters. These parameters are different depending on the method selected. 
#' Every parameter has a default value set in the atomic functions in case that no selection is wanted. 
#' Everything concerning these parameters is explained in the section \code{Details} of the function \code{\link[downscaleR]{downscaleTrain}}. However, if wanted, the atomic functions can be seen here: 
#' \code{\link[downscaleR]{glm.train}} and \code{\link[deepnet]{nn.train}}.  
#' @param path A string indicating the path where to save the prediction. 
#' @return Saves the prediction where specified.
#' @seealso 
#' downscaleTrain for training a downscaling model
#' downscalePredict for prediction for a a test dataset with a trained model for 
#' \href{https://github.com/SantanderMetGroup/downscaleR/wiki/training-downscaling-models}{downscaleR Wiki} for downscaling seasonal forecasting and climate projections.
#' @family downscaling.functions
#' @author J. Bano-Medina
#' @export

downscaleChunk <- function(x, y, newdata, simulate = FALSE,
                           method, ...,
                           prepareData.args = list("global.vars" = NULL, "combined.only" = TRUE, "spatial.predictors" = NULL, "local.predictors" = NULL, "extended.predictors" = NULL),
                           condition = NULL, threshold = NULL, predict = TRUE,
                           path = getwd()) {
  
  if (!exists("global.vars",prepareData.args)) prepareData.args$global.vars <- NULL
  if (!exists("combined.only",prepareData.args)) prepareData.args$combined.only <- TRUE
  if (!exists("spatial.predictors",prepareData.args)) prepareData.args$spatial.predictors <- NULL
  if (!exists("local.predictors",prepareData.args)) prepareData.args$local.predictors <- NULL
  if (!exists("extended.predictors",prepareData.args)) prepareData.args$extended.predictors <- NULL
  
  x <- getTemporalIntersection(x,y,which.return = "obs")
  y <- getTemporalIntersection(x,y,which.return = "prd")
  
  latitudes <- getCoordinates(y)$y
  chunks <- length(latitudes)
  lapply(1:chunks,FUN = function(z){
    print(paste("Training chunk:",z,"out of",chunks))
    y_chunk <- subsetDimension(y,dimension = "lat", indices = z)
    xyT <- prepareData(x = x, y = y_chunk, global.vars = prepareData.args$global.vars, combined.only = prepareData.args$combined.only, spatial.predictors = prepareData.args$spatial.predictors, local.predictors = prepareData.args$local.predictors, extended.predictors = prepareData.args$extended.predictors)
    model <- downscaleTrain(xyT, method = method, condition = condition, threshold = threshold, predict = predict, model.verbose = FALSE, ...)
    
    p <- lapply(simulate, function(sim) {
      lapply(newdata, function(zz) {
        xyt <- prepareNewData(zz,xyT)
        downscalePredict(newdata = xyt,model, simulate = sim)
      })
    }) %>% unlist(recursive = FALSE)
    
    if (z < 10) {zn <- paste0("00",z)}
    else if (z < 100 & z >= 10) {zn <- paste0("0",z)}
    else{zn <- z}
    ini <- ifelse(isTRUE(predict),length(p)+1,length(p))
    lapply(1:ini, function(zzz) {
      if (ini == (length(p)+1) & zzz == 1) {
        grid <- model$pred
      } else if (ini == (length(p)+1) & zzz != 1) {
        grid <- p[[zzz-1]] 
      } else if (ini == (length(p))) {
        grid <- p[[zzz]]
      }
      save(grid, file = paste0(path,"/dataset",zzz,"_chunk",zn,".rda"))
    })
    p <- NULL
  })
  ini <- ifelse(isTRUE(predict),length(newdata)+1,length(newdata))
  pred <- list()
  for (i in 1:ini) {
    lf <- list.files("./", pattern = paste0("dataset",i), full.names = TRUE)
    chunk.list <- lapply(lf, function(x) get(load(x)))
    pred[[i]] <- bindGrid(chunk.list, dimension = "lat")
    file.remove(lf)
  }
  if (ini == 1) pred <- pred[[1]]
  return(pred)
}
