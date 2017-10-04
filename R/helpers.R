#' @title Compute dates of a downscaled observational dataset
#' @description The function calculates the appropiate Dates slot in the returned output of downscaling functions,
#' considering the possible mismatches in time resolution between predictors and predictand, the multigrid dates slot etc.
#' @param obs.dates The \code{Dates} slot of the 'obs' input in the downscaling method
#' @param sim.dates The \code{Dates} slot of the 'sim' input in the downscaling method
#' @return A new \code{Dates} list that preserves the temporal extent of the downscaled simulations but considering
#' the temporal resolution of the downscaled variable
#' @details The function is intended for internal use only. Sometimes the time resolution of the predictors does not match 
#' that of the downscaled variable (e.g., suppose that instantaneous surface temperature at 12:00 UTC is used as predictor
#'  of daily minimum temperature). In addition, in case of multiple predictors the \code{Dates} slot of the simulated series
#'  has several start/end time lists, one for each predictor, while there is only one predictand. For this reason,
#'  the function takes care of adjusting adequately the returned \code{Dates} slot.
#'  @author J. Bedia 
#'  @keywords internal
#'  @export

dateReplacement <- function(obs.dates, sim.dates) {
    time.res <- difftime(as.POSIXlt(obs.dates$end[1]), as.POSIXlt(obs.dates$start[1]))
    hours <- as.POSIXlt(obs.dates$start)$hour
    tz <- tryCatch(unlist(strsplit(obs.dates$start[1], split = "\\s"))[3], error = function(er) {er <- ""})
    if (is.null(names(sim.dates))) {
        tz.sim <- unlist(strsplit(sim.dates[[1]]$start[1], split = "\\s"))[3]
        sim.dates.ref <- as.POSIXlt(sim.dates[[1]]$start, tz = tz.sim)
    } else {
        tz.sim <- unlist(strsplit(sim.dates$start[1], split = "\\s"))[3]
        sim.dates.ref <- as.POSIXlt(sim.dates$start, tz = tz.sim)
    }
    aux.string <- paste(sim.dates.ref$year + 1900, sim.dates.ref$mon + 1, sim.dates.ref$mday, hours, sep = "-")
    length(aux.string) <- length(sim.dates.ref)
    start <- strptime(aux.string, "%Y-%m-%d-%H", tz)
    aux.string <- NULL
    end <- as.POSIXct(start + time.res)
    start <- as.POSIXct(start)
    usetz <- ifelse(identical(tz, ""), FALSE, TRUE)
    start <- format.POSIXct(start, "%Y-%m-%d %H:%M:%S", usetz = usetz)
    end <- format.POSIXct(end, "%Y-%m-%d %H:%M:%S", usetz = usetz)
    return(list("start" = start, "end" = end))
}
# End



#' @title Set the 'dimensions' attribute 
#' @description Sets the 'dimensions' attribute of model out Data objects after downscaling
#' @param obs A observations object
#' @param multi.member Logical indicating if simulation data is a multimember
#' @return A character vector indicating the dimensions of the output object
#' @keywords internal
#' @importFrom transformeR getDim
#' @author J. Bedia


renameDims <- function(obs, multi.member) {
    dimNames <- getDim(obs)
    # Remove "station" from dimensions for single-station objects
    st.dim.index <- grep("loc", dimNames)
    if (!identical(st.dim.index, integer(0))) {
        dim.st <- dim(obs$Data)[st.dim.index]
        if (identical(dim.st, 1L)) {
            dimNames <- dimNames[-st.dim.index]
        }
    }
    if (isTRUE(multi.member)) dimNames <- c("member", dimNames)
    return(dimNames)
}
# End


#' @title getIntersect
#' @description Get the common period of the objects obs and prd
#' @author S. Herrera
#' @keywords internal
#' @importFrom transformeR subsetDimension

getIntersect <- function(obs,prd){
    dimNames <- attr(obs$Data, "dimensions")
    indDates <- which(as.POSIXct(obs$Dates$start, tz = "GMT", format = "%Y-%m-%d") == as.POSIXct(prd$Dates$start, tz = "GMT", format = "%Y-%m-%d"))
    auxDates <- as.POSIXct(obs$Dates$start[indDates], tz = "GMT", format = "%Y-%m-%d")
    indObs <- which(is.element(as.POSIXct(obs$Dates$start, tz = "GMT", format = "%Y-%m-%d"), auxDates))
    obs <- subsetDimension(obs, dimension = "time", indices = indObs)
    dimNames <- attr(prd$Data, "dimensions")
    indObs <- which(is.element(as.POSIXct(prd$Dates$start, tz = "GMT", format = "%Y-%m-%d"), auxDates))
    prd <- subsetDimension(prd, dimension = "time", indices = indObs)
    obj <- list(obs = obs, prd = prd)
    obj$Dates$start <- as.POSIXct(obs$Dates$start, tz = "GMT", format = "%Y-%m-%d")
    obj$Dates$end <- as.POSIXct(obs$Dates$end, tz = "GMT", format = "%Y-%m-%d")
    attr(obj$obs$Data, "dimensions") <- attr(obs$Data, "dimensions")
    attr(obj$prd$Data, "dimensions") <- attr(prd$Data, "dimensions")
    return(obj)
}


#' @title Get grid coordinates as 2D matrix
#' @description Obtain grid coordinates as 2D matrix
#' @param grid An input grid
#' @return A 2D matrix of x-y coordinates (in this order)
#' @importFrom transformeR typeofGrid
#' @keywords internal
#' @author J. Bedia

get2DmatCoordinates <- function(grid) {
    if (typeofGrid(grid) == "regular_grid") {
        getCoordinates(grid) %>%  expand.grid()
    } else if (typeofGrid(grid) == "station") {
        getCoordinates(grid)
    } else if (typeofGrid(grid) == "rotated_grid") {
        stop("Direct downscaling of rotated grids is not supported", call. = FALSE)
    }
}


#' @title Obtain a grob object from ordinary plot
#' @description Obtain a grob object from ordinary plot
#' @return A grob object
#' @keywords internal
#' @importFrom gridGraphics grid.echo 
#' @importFrom grid grid.grab
#' @author J. Bedia


grabGrob <- function(){
      grid.echo()
      grid.grab()
}
# End