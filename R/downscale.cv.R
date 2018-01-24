##############################################################################################################
#                     GENERAL DOWNSCALING                                                                    #
##############################################################################################################
##     downscale.cv.R Downscale climate with cross validation.
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

#' @title Downscale climate data and reconstruct the temporal serie by splitting the data in k folds.
#' @description Downscale climate data and reconstruct the temporal serie by splitting the data in k folds. 
#' Statistical downscaling methods are: analogs, generalized linear models (GLM) and Neural Networks (NN). 
#' @param x The input grid. It should be an object as returned by \pkg{loadeR}.
#' @param y The observations dataset. It should be an object as returned by \pkg{loadeR}.
#' @param method A string value. Type of transer function. Options are c("analogs","GLM","NN").
#' @param folds Could be a fraction, value between (0,1) indicating the fraction of the data that will define the train set, 
#' or an integer indicating the number of folds. It can also be a list of folds indicating the years of each fold. 
#' @param type A string, c("chronological","random"). Indicates how to split the data in folds. Default is "chronological".
#' @param scale A logical value. If TRUE then the data is standardized. Default is FALSE.
#' @param singlesite A logical value. Wether to perform the study singlesite or multisite. Multisite option is only available when
#' the selected method is or analogs or NN. For GLM, multisite can only be performed when the optional parameter of GLM's \code{fitting}, is fitting = "MP". 
#' @param global.vars An optional character vector with the short names of the variables of the input x multigrid to be retained as global predictors 
#' (use the getVarNames helper if not sure about variable names). 
#' This argument just produces a call to subsetGrid, but it is included here for better flexibility in downscaling experiments (predictor screening...). 
#' For instance, it allows to use some specific variables contained in x as local predictors and the remaining ones, specified in subset.vars, 
#' as either raw global predictors or to construct the combined PC.
#' @param PCA Default to NULL, and not used. Otherwise, a named list of arguments in the form argument = value, 
#' with the arguments to be passed to prinComp to perform Principal Component Analysis of the predictors grid (x). 
#' See Details on principal component analysis of predictors.
#' @param combined.only Optional, and only used if PCA parameters are passed. Should the combined PC be used as the only global predictor? Default to TRUE. 
#' Otherwise, the combined PC constructed with which.combine argument in prinComp is append to the PCs of the remaining variables within the grid.
#' @param local.predictors Default to \code{NULL}, and not used. Otherwise, a named list of arguments in the form \code{argument = value},
#'  with the following arguments:
#'  \itemize{
#'    \item \code{neigh.vars}: names of the variables in \code{x} to be used as local predictors
#'    \item \code{neigh.fun}: Optional. Aggregation function for the selected local neighbours.
#'    The aggregation function is specified as a list, indicating the name of the aggregation function in
#'     first place (as character), and other optional arguments to be passed to the aggregation function.
#'     For instance, to compute the average skipping missing values: \code{neigh.fun = list(FUN= "mean", na.rm = TRUE)}.
#'     Default to NULL, meaning that no aggregation is performed.
#'    \item \code{n.neighs}: Integer. Number of nearest neighbours to use. If a single value is introduced, and there is more
#'    than one variable in \code{neigh.vars}, the same value is used for all variables. Otherwise, this should be a vector of the same
#'    length as \code{neigh.vars} to indicate a different number of nearest neighbours for different variables.
#'  }
#' @param neurons A numeric value. Indicates the size of the random nonlinear dimension where the input data is projected.
#' @param filt A logical expression (i.e. = ">0"). This will filter all values that do not accomplish that logical statement. Default is NULL.
#' @param ... Optional parameters. These parameters are different depending on the method selected. 
#' Every parameter has a default value set in the atomic functions in case that no selection is wanted. 
#' Everything concerning these parameters is explained in the section \code{Details} of the function \code{\link[downscaleR]{downscale.train}}. However, if wanted, the atomic functions can be seen here: 
#' \code{\link[downscaleR]{glm.train}} and \code{\link[deepnet]{nn.train}}.  
#' @details The functon relies on \code{\link[downscaleR]{prepare_predictors}}, \code{\link[downscaleR]{prepare_newdata}}, \code{\link[downscaleR]{downscale.train}}, and \code{\link[downscaleR]{downscale.predict}}. 
#' For more information please visit these functions.
#' If the variable to downscale is the precipitation and it is a binary variable, then two temporal series will be returned:
#' 1) The temporal serie with binary values filtered by a threshold adjusted by the train dataset, see \code{\link[transformeR]{convert2bin}} for more details.
#' 2) The temporal serie with the results obtained by the downscaling, without any binary converting process.
#' We recommend to get rid of the NaN/NA when dealing with multisite mode.
#' 
#' @return The reconstructed downscaled temporal serie.
#' @importFrom transformeR dataSplit convert2bin localScaling
#' @author J. Bano-Medina
#' @export
#' @examples 
#' x <- makeMultiGrid(NCEP_Iberia_hus850, NCEP_Iberia_ta850)
#' x <- subsetGrid(x, years = 1985:1995)
#' # Loading predictands
#' y <- VALUE_Iberia_pr
#' y <- getTemporalIntersection(obs = y,prd = x, "obs" )
#' x <- getTemporalIntersection(obs = y,prd = x, "prd" )
#' # Reconstructing the downscaled serie in 3 folds
#' pred <- downscale.cv(x,y,folds = 3,type = "chronological", 
#'         scale = TRUE, method = "GLM", filt = ">0")
#' # ... or with dates ...
#' pred <- downscale.cv(x,y,type = "chronological", scale = TRUE, 
#'                      method = "GLM", filt = ">0",
#'                      folds = list(c("1985","1986","1987","1988"),
#'                                   c("1989","1990","1991","1992"),
#'                                   c("1993","1994","1995")))
#' # Reconstructing the downscaled serie in 3 folds with a 
#' # pre-processed of the predictors with principal component analysis.
#' pred <- downscale.cv(x,y,folds = 3,type = "chronological", 
#'                       method = "GLM", family = Gamma(link = "log"), filt = ">0",
#'                      PCA = list(which.combine = getVarNames(x),v.exp = 0.9))
#' # Reconstructing the downscaled serie in 3 folds with local predictors.
#' pred <- downscale.cv(x,y,folds = 3,type = "chronological", 
#'                      scale = TRUE, method = "GLM", filt = ">0",
#'                      local.predictors = list(neigh.vars = "shum@850",n.neighs = 4))

downscale.cv <- function(x, y, method,
                         folds = 4, type = "chronological", scale = FALSE, singlesite = TRUE,
                         global.vars = NULL, PCA = NULL, combined.only = TRUE, local.predictors = NULL, neurons = NULL,
                         filt = NULL, ...) {
  data <- dataSplit(x,y, f = folds, type = type)
  p <- lapply(1:length(data), FUN = function(xx) {
    print(paste("fold:",xx,"-->","calculating..."))
    xT <- data[[xx]]$train$x ; yT <- data[[xx]]$train$y
    xt <- data[[xx]]$test$x  ; yt <- data[[xx]]$test$y
    if (isTRUE(scale)) {
      xt <- localScaling(xt, base = xT, scale = scale)
      xT <- localScaling(xT, base = xT, scale = scale)
    }
    xT <- prepare_predictors(x = xT, y = yT, global.vars, PCA, combined.only, local.predictors, neurons)
    xt <- prepare_newdata(newdata = xt, predictor = xT)
    model <- downscale.train(xT, method, singlesite, filt, ...)
    if (all(as.vector(y$Data) %in% c(0,1,NA,NaN), na.rm = TRUE)) {
      y.prob <- downscale.predict(xt, model)
      if (singlesite == TRUE && method == "GLM") {
        if (model$conf$atomic_model[[1]]$info$simulate == "yes") {
          y.bin  <- y.prob}
        else {
          y.bin  <- convert2bin(y.prob, ref.obs = yT, ref.pred = model$pred)}}
      else{
        y.bin  <- convert2bin(y.prob, ref.obs = yT, ref.pred = model$pred)}
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
