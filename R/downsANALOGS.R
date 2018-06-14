##############################################################################################################
#                     Analogs                                                                                #
##############################################################################################################
#' @title Analogs
#' @description Analog method implementation
#' @param x The grid data. Class: matrix.
#' @param y The observations data. Class: matrix.
#' @param dates Dates of the grid and observations data.
#' @param n.analogs An integer. Number of analogs. Default is 4.
#' @param sel.fun A string. Select a function to apply to the analogs selected for a given observation. Options are 
#' "mean", "wmean" (i.e., weighted mean), "max", "min", "median", "prcXX" 
#' (i.e., prc85 means the 85th percentile of the analogs values distribution). Default is "mean".
#' the function applied to the analogs values, (i.e., sel.fun = c("mean","max","min","median","prcXX"), with default "mean") 
#' and the temporal window, (i.e., window = 0).
#' @param window An integer. Window of days removed when selecting analogs. 
#' If window = 7, then 7 days after the observation date and the 7 days before the observation date are removed. Default is 0.
#' @param n.random An integer. Choose N random analogs among the closest n.analogs. Default is NULL.
#' @param pool An integer. Number of auxiliary analogs in case there are NaN or NA in the original analogs.
#' @return A list containing the grid data, the observations and a list with information concerning the analogs.
#' @details The analogs actually is not a model in the sense that it is not really trained. 
#' The model would be the data where to search the analogs. For this reason this function only saves the information to search analogs in a list. 
#' @author J. Bano-Medina
analogs.train <- function(x, y, dates, n.analogs = 4, sel.fun = "mean", window = 7, n.random = NULL, pool = 0){
  if (n.analogs == 1) {sel.fun <- NULL}
  return(list("dataset_x" = x, "dataset_y" = y, "dates" = list("train" = dates), "info" = list("n.analogs" = n.analogs, "sel.fun" = sel.fun, "window" = window, "n.random" = n.random, "pool" = pool)))}

#' @title Donwscaling with the analogs method (test).
#' @description Donwscaling with the analogs method (test).
#' @param newdata The grid data.  Class: matrix.
#' @param x The grid data where to search the analogs. Class: matrix.
#' @param y The observations data of the same days of x. Class: matrix.
#' @param dates A list containing both the newdata grid  and x grid, dates.
#' @param info A list containing the information required to perform the analogs method: number of analogs, function, window. These parameters are grouped in a list in \code{\link[downscaleR]{analogs.train}}.
#' @return A matrix containing the predictions.
#' @details The selected functions use the base R functions: \code{\link[base]{mean}}, \code{\link[stats]{median}}, \code{\link[stats]{weighted.mean}}, \code{\link[stats]{quantile}}, \code{\link[base]{Extremes}}.
#' @author J. Bano-Medina
#' @importFrom stats dist weighted.mean
#' @importFrom fields rdist
analogs.test <- function(newdata, x, y, dates, info) {
  # Dealing with dates
  dist2test <- rdist(x,newdata) 
  date_newdata <- julian(as.Date(dates$test),origin = as.Date("1970-01-01"))
  date_x <- julian(as.Date(dates$train),origin = as.Date("1970-01-01"))
  zeros_window <- apply(as.matrix(date_newdata),MARGIN = 1,function(x){
    date_window <- seq(from = x - info$window, to = x + info$window, by = 1)
    match(date_window,date_x)
  })
  num_zeros <- apply(as.matrix(date_newdata),MARGIN = 1,function(x){
    date_window <- seq(from = x - info$window, to = x + info$window, by = 1)
    ind_window <- match(date_window,date_x)
    length(which(is.na(ind_window) == FALSE)) 
  })
  # Searching analogs
  zeros_window <- matrix(zeros_window,nrow = length(-info$window:info$window),ncol = ncol(dist2test))
  num_zeros <- matrix(num_zeros,nrow = 1,ncol = ncol(dist2test))
  apply(rbind(dist2test,zeros_window,num_zeros),MARGIN = 2, FUN = function(x){
    ind_zeros <- x[length(x)]
    x[x[((length(x) - length(-info$window:info$window)):(length(x) - 1))]] <- 0
    x <- x[-((length(x) - length(-info$window:info$window)):length(x))]
    dist.sorted <- sort(x)
    if (ind_zeros != 0) {dist.analogs <- matrix(rep(dist.sorted[-(1:ind_zeros)][1:info$n.analogs],ncol(y)),nrow = info$n.analogs,ncol = ncol(y))}
    else{dist.analogs <- matrix(rep(dist.sorted[1:info$n.analogs],ncol(y)),nrow = info$n.analogs,ncol = ncol(y))}
    ind.analogs <- unique(match(dist.analogs,x))
    value.analogs <- matrix(y[ind.analogs,], nrow = info$n.analogs, ncol = ncol(y))
    # Pool analogs (if wanted and if necessary..)
    if (anyNA(value.analogs) && (info$pool != 0)) {
      if (ind_zeros != 0) {dist.pool <- dist.sorted[-(1:ind_zeros)][(info$n.analogs + 1):(info$n.analogs + info$pool)]}
      else{dist.pool <- dist.sorted[(info$n.analogs + 1):(info$n.analogs + info$pool)]}
      ind.pool <- match(dist.pool,x)
      value.pool <- matrix(y[ind.pool,], nrow = info$pool, ncol = ncol(y))
      value.analogs <- sapply(1:ncol(value.analogs), function(i) {
        if (anyNA(value.analogs[,i])) {
          ind.no <- which(is.na(value.analogs[,i]))
          ind.yes <- setdiff(1:nrow(value.analogs),ind.no)
          value.analogs[,i] <- c(value.analogs[ind.yes,i],value.analogs[ind.no,i])
          limit_replacement2 <- min(c(length(ind.no),nrow(value.pool)))
          limit_replacement1 <- nrow(value.analogs) - length(ind.no) + limit_replacement2
          value.analogs[(nrow(value.analogs) + 1 - length(ind.no)):limit_replacement1,i] <- value.pool[1:limit_replacement2,i]
          value.analogs[,i]}
        else{
          value.analogs[,i]}
      })
      value.analogs <- matrix(value.analogs, nrow = info$n.analogs, ncol = ncol(y))
      dist.analogs <- sapply(1:ncol(value.analogs), function(i) {
        if (anyNA(value.analogs[,i])) {
          ind.no <- which(is.na(value.analogs[,i]))
          ind.yes <- setdiff(1:nrow(value.analogs),ind.no)
          value.analogs[,i] <- c(value.analogs[ind.yes,i],value.analogs[ind.no,i])
          limit_replacement2 <- min(c(length(ind.no),nrow(value.pool)))
          limit_replacement1 <- nrow(value.analogs) - length(ind.no) + limit_replacement2
          dist.analogs[(nrow(value.analogs) + 1 - length(ind.no)):limit_replacement1,i] <- dist.pool[1:limit_replacement2]
          dist.analogs[,i]}
        else{
          dist.analogs[,i]}
      })
    }
    dist.analogs <- matrix(dist.analogs, nrow = info$n.analogs, ncol = ncol(y))
    # Random analogs (if selected)
    if (!is.null(info$n.random)) {
      ind.random <- sample(1:dim(value.analogs)[1],size = info$n.random, replace = FALSE)
      dist.analogs  <- dist.analogs[ind.random,]
      value.analogs <- value.analogs[ind.random,]}
    # Selection function applied to the ensemble of analogs
    if (is.null(info$sel.fun)) {
      value.analogs}
    else if (info$sel.fun == "mean" || info$sel.fun == "median" || info$sel.fun == "max" || info$sel.fun == "min") {
      apply(X = value.analogs, MARGIN = 2, FUN = function(X){do.call(args = list(X,"na.rm" = TRUE), what = info$sel.fun)})}
    else if (substr(info$sel.fun, 1, 3) == "prc") {
      percentile <- as.numeric(substr(info$sel.fun, 4, 5))/100
      apply(X = value.analogs, MARGIN = 2, FUN = function(X){do.call(args = list(X,percentile,"na.rm" = TRUE), what = "quantile")})}
    else if (info$sel.fun == "wmean") {
      w <- 1/dist.analogs
      w <- 1/dist.analogs
      sapply(1:ncol(y), function(zz){
      weighted.mean(value.analogs[,zz],w[,zz],na.rm = TRUE)})}
    else {
      warning("Unknown selected function")}
  }) %>% t() %>% drop() 
}


