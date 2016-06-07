#' @title Parallel availability check
#' @description Check availability of parallel package and set parameters
#' @param parallel Logical, should parallel execution be used?
#' @param maxcores upper bound for self-selected number of cores
#' @param ncores number of cores used in parallel computation, self-selected number of cores
#'  is used when \code{is.null(ncpus)} (the default).
#' @details The function checks if package parallel is available. Then, it is checked whether the FORK nodes
#' can be initialized
#' @return A list with two elements:
#' \itemize{
#' \item \code{hasparallel}, a logical flag indicating if parallelization is enabled and
#' \item \code{cl}: parallel socket cluster object, or NULL.
#' }
#' @author J. Bedia, relying on previous code by Jonas Bhend and Matteo de Felice
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @keywords internal

parallelCheck <- function(parallel, max.ncores = 16, ncores = NULL) {
      hasparallel <- FALSE
      .cl <- NULL
      if (parallel & grepl("windows", .Platform$OS.type, ignore.case = TRUE)) {
            message("Parallelization is not supported on Windows machines")    
      } 
      if (parallel && requireNamespace("parallel", quietly = TRUE)) {
            if (is.null(ncores)) {
                  ncores <- min(max(parallel::detectCores() - 1, 1), max.ncores)
            }
            if (ncores > 1) {
                  .cl <- try(parallel::makeCluster(ncores, type = 'FORK'), silent = TRUE)
                  if (!"try-error" %in% class(.cl)) {
                        hasparallel <- TRUE
                        message("Parallel computing enabled\nNumber of workers: ", ncores)
                  } else {
                        .cl <- NULL
                  }
            }
      } else if (parallel) {
            message("Parallel computing not enabled (parallel package is not available)")
      }
      list("hasparallel" = hasparallel, "cl" = .cl)
}
