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
  df <- data.frame(cbind(y,x)); colnames(df) <- paste0('X',1:(dim(x)[2] + 1))
  weights <- glm(X1~.,data = df, family)
  return(weights)}

#' @title Donwscaling with a given generalized linear model (GLM).
#' @description Donwscaling with a given generalized linear models (GLM) calculated in \code{\link[downscaleR]{downscale.predict}} via \code{\link[downscaleR]{glm.train}}.
#' @param x The grid data. Class: matrix.
#' @param weights Object as returned from \code{\link[downscaleR]{glm.train}}
#' @return The predicted matrix.
#' @details Predicts by using the base function \code{\link[stats]{predict}}. This function is internal and should not be used by the user. The user should use \code{\link[downscaleR]{downscale.predict}}.
#' @author J. Bano-Medina
#' @export
glm.predict <- function(x, weights){
  df <- data.frame(x); colnames(df) <- paste0('X',2:(dim(x)[2] + 1))
  if (weights$family$family == "binomial") {
    pred <- predict(weights, newdata = df, type = 'response')}
  else{
    pred <- predict(weights, newdata = df)}
  return(pred)}