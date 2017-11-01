##############################################################################################################
#                     Analogs                                                                                #
##############################################################################################################
#' @title Donwscaling with the analogs method (train).
#' @description Donwscaling with the analogs method (train).
#' @param x The grid data. Class: matrix.
#' @param y The observations data. Class: matrix.
#' @param n.analogs Number of analogs. Default is 4.
#' @param sel.fun Function applied to the analogs. Can be "mean", "wmean" (i.e., weighted mean), "max", "min", "median", "prcXX" (i.e., prc85 means the 85th percentile of the analogs values distribution). Default is "mean".
#' @param window Window of days removed when selecting analogs. Default is 7. This means that the 7 days after the date and the 7 days before the date are removed.
#' @return A list containing the grid data, the observations and a list with information concerning the analogs.
#' @details The analogs actually is not a model in the sense that it is not really trained. The model would be the data where to search the analogs. For this reason this function only saves the information to search analogs in a list. 
#' @author J. Bano-Medina
#' @export
analogs.train <- function(x, y, n.analogs = 4, sel.fun = "mean", window = 0){
  #if (window == 0) {warning("Window is 0, this may lead to an unrealistic prediction.")}
  if (n.analogs == 1) {sel.fun <- NULL}
  return(list("dataset_x" = x, "dataset_y" = y, "info" = list("n.analogs" = n.analogs, "sel.fun" = sel.fun, "window" = window)))}

#' @title Donwscaling with the analogs method (test).
#' @description Donwscaling with the analogs method (test).
#' @param newdata The grid data.  Class: matrix.
#' @param x The grid data where to search the analogs. Class: matrix.
#' @param y The observations data of the same days of x. Class: matrix.
#' @param info A list containing the information required to perform the analogs method: number of analogs, function, window. These parameters are grouped in a list in \code{\link[downscaleR]{analogs.train}}.
#' @return A matrix containing the predictions.
#' @details The selected functions use the base R functions: \code{\link[base]{mean}}, \code{\link[stats]{median}}, \code{\link[stats]{weighted.mean}}, \code{\link[stats]{quantile}}, \code{\link[base]{Extremes}}.
#' @author J. Bano-Medina
#' @importFrom stats dist
#' @export
analogs.test <- function(newdata, x, y, info){
  prediction <- sapply(1:dim(newdata)[1], FUN = function(z){
    print(z)
    if (info$window == 0) {
      ind <- z}
    else {
      if (z <= info$window) {
        ind <- 1:(z + info$window)}
      else if (z + info$window > dim(x)[1]) {
        ind <- (z - info$window):(z + dim(x)[1])}
      else {
        ind <- (z - info$window):(z + info$window)}}
    ind_yes <- setdiff(1:dim(newdata)[1],ind)
    dist2test <- sapply(1:dim(newdata)[1], FUN = function(xx) {
      if (sum(xx == ind) == 1) {0}
      else {dist(rbind(newdata[z,],x[ind_yes[xx],]))}})
    ind_zeros <- which(dist2test == 0)
    dist.analogs <- sort(dist2test[-ind_zeros])[1:info$n.analogs]
    ind.analogs <- sapply(1:info$n.analogs, FUN = function(zz){
      ind.analogs <- which((dist2test == dist.analogs[zz]) == TRUE)})
    value.analogs <- y[ind.analogs,]
    if (is.null(info$sel.fun)) {
      pred <- value.analogs}
    else if (info$sel.fun == "mean" || info$sel.fun == "median" || info$sel.fun == "max" || info$sel.fun == "min") {
      pred <- apply(X = value.analogs, MARGIN = 2, FUN = function(X){do.call(args = list(X), what = info$sel.fun)})}
    else if (substr(info$sel.fun, 1, 3) == "prc") {
      percentile <- as.numeric(substr(info$sel.fun, 4, 5))/100
      pred <- apply(X = value.analogs, MARGIN = 2, FUN = function(X){do.call(args = list(X,percentile), what = "quantile")})}
    else if (info$sel.fun == "wmean") {
      w <- 1/dist.analogs
      pred <- apply(X = value.analogs, MARGIN = 2, FUN = function(X){do.call(args = list(X,w), what = "weighted.mean")})}
    else {
      warning("Unknown selected function")}
    return(pred)})
  return(t(prediction))
}


