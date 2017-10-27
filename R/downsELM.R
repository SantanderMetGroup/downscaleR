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