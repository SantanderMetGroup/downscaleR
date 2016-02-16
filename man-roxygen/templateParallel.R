#' @section Parallel Processing:
#' 
#' Parallel processing is enabled using the \pkg{parallel} package. 
#' Parallelization is undertaken by a FORK-type parallel socket cluster formed by \code{ncores}.
#' If \code{ncores} is not specified (default), \code{ncores} will be one less than the autodetected number of cores.
#' The maximum number of cores used for parallel processing can be set with the \code{max.ncores} argument, 
#' although this will be reset to the auto-detected number of cores minus 1 if this number is exceeded. Note that not all 
#' code, but just some critical loops within the function are parallelized.
#' 
#' In practice, parallelization does not always result in smaller execution times, due to the parallel overhead.
#' However, parallel computing may potentially provide a significant speedup for the 
#' particular case of large multimember datasets or large grids.
#'  
#' Parallel computing is currently not available for Windows machines.

