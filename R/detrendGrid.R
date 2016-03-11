#' @title Linear detrending
#' @description Perform a linear detrending along the time dimension of a grid
#' @param grid Input grid (possibly multimember)
#' @template templateParallelParams
#' @return A detrended grid. Values are residuals.
#' @details  Performs a simple linear detrending by fitting a linear model and retaining the residuals.
#' An attribute indicating the linear detrending is added to the \code{Variable} component of the output grid.
#' @template templateParallel
#' @export
#' @importFrom parallel stopCluster
#' @importFrom parallel parApply
#' @export
#' @author J Bedia

detrendGrid <- function(grid, parallel = FALSE, max.ncores = 16, ncores = NULL) {
      arr <- grid$Data
      refdim <- dim(arr)
      dimNames <- attr(arr, "dimensions")
      # grid$Data <- NULL
      mar <- grep("time", attr(arr, "dimensions"), invert = TRUE)
      ntimes <- dim(arr)[grep("time", attr(arr, "dimensions"))]
      x <- 1:ntimes
      message("[", Sys.time(), "] - Detrending...")
      arr <- if (isTRUE(parallel)) {
            parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
            on.exit(parallel::stopCluster(parallel.pars$cl))
            unname(parallel::parApply(cl = parallel.pars$cl, arr, MARGIN = mar, FUN = function(y) {
                  out <- rep(NA, length(y))
                  ind <- intersect(which(!is.na(y)), which(!is.na(x)))
                  out[ind] <- tryCatch(expr = summary(lm(y ~ I(x)))$resid,
                                       error = function(err) {
                                             rep(NA,ntimes)
                                       })
            })
            )
      } else {
            unname(apply(arr, MARGIN = mar, FUN = function(y) {
                  out <- rep(NA, length(y))
                  ind <- intersect(which(!is.na(y)), which(!is.na(x)))
                  out[ind] <- tryCatch(expr = summary(lm(y ~ I(x)))$resid,
                                       error = function(err) {
                                             rep(NA,ntimes)
                                       })
            })
            )
      }
      message("[", Sys.time(), "] - Done.")
      newdim <- dim(arr)
      if (any(newdim != refdim)) arr <- aperm(arr, perm = match(newdim, refdim))
      grid$Data <- arr
      attr(grid$Data, "dimensions") <- dimNames
      attr(grid$Variable, "detrended:method") <- "linear"
      attr(grid$Variable, "units")
      return(grid)
}
# End
