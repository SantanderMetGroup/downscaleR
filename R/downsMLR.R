##############################################################################################################
#                     Multiple Linear Regression (GLM)                                                       #
##############################################################################################################
#' @title Donwscaling with multiple linear regression (MLR).
#' @description Donwscaling with multiple linear regression (MLR).
#' @param x The grid data. Class: matrix.
#' @param y The observations data. Class: matrix.
#' @param fitting  Type of fitting used: c("LS","MP","MP+L2"). There are three options: least squares ("LS") that downscales via \code{\link[stats]{lm}}, Moore-Penrose ("MP") that uses the generalized inverse via \code{\link[MASS]{ginv}} or the L2 penalty, commonly known as ridge regression ("MP+L2"). This is an optional parameter and the default is "LS". 
#' @return If fitting = "LS", then returns a MLR model as returned from \code{\link[stats]{lm}} whereas if fitting = ("MP" or "MP+L2") then returns just the estimated coefficients/weights.
#' @details If "MP" or "MP+L2" are selected, the bias is input as an additional input neuron in the weight matrix with value equal 1. The indicated optional parameters can be changed by including them in \code{\link[downscaleR]{downscale.predict}}, in this case \code{fitting} is the only optional parameter. This function is internal and should not be used by the user. The user should use \code{\link[downscaleR]{downscale.predict}}.
#' @author J. Bano-Medina
#' @export
mlr.train <- function(x, y, fitting = 'LS'){
  if (fitting == 'LS') {
    df <- data.frame(cbind(y,x)); colnames(df) <- paste0('X',1:(dim(x)[2] + 1))
    weights <- lm(X1~.,data = df)}
  else{
    neurons <- vector("integer",2)
    neurons[1] <- dim(x)[2] + 1
    if (is.null(dim(y)[2])) {neurons[2] <- 1}
    else{neurons[2] <- dim(y)[2]}
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
      weights <- solve((t(data.matrix(x)) %*% data.matrix(x) + (diag(dim(data.matrix(x))[2]) / C_winner))) %*% (t(data.matrix(x)) %*% data.matrix(y))}}
  return(list("weights" = weights, "info" = list("fitting" = fitting)))}

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
    df <- data.frame(x); colnames(df) <- paste0('X',2:(dim(x)[2] + 1))
    pred <- predict(weights, newdata = df, type = 'response')}
  if (info$fitting == 'MP' || info$fitting == 'MP+L2') {
    pred <- data.matrix(x) %*% weights}
  return(pred)}
