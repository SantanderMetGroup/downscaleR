#' @title Compute dates of a downscaled observational dataset
#' @description The function calculates the appropiate Dates slot in the returned output of downscaling functions,
#' considering the possible mismatches between predictors and predictand, the multifield dates slot etc.
#' @param obs.dates The \code{Dates} slot of the 'obs' input in the downscaling method
#' @param sim.dates The \code{Dates} slot of the 'sim' input in the downscaling method
#' @return A new \code{Dates} list that preserves the temporal extent of the downscaled simulations but considering
#' the temporal resolution of the downscaled variable
#' @details The function is intended for internal use only. Sometimes the time resolution of the predictors does not match 
#' that of the downscaled variable (e.g., suppose that instantaneous surface temperature at 12:00 UTC is used as predictor
#'  of daily minimum temperature). In addition, in case of multiple predictors the \code{Dates} slot of the simulated series
#'  has several start/end time lists, one for each predictor, while there is only one predictand. For this reason,
#'  the function takes care of adjusting adequately the returned \code{Dates} slot.
#'  @author J. Bedia \email{joaquin.bedia@@gmail.com}
#'  @keywords internal
#'  @export
#' 

dateReplacement <- function(obs.dates, sim.dates) { 
      time.step <- diff(as.POSIXlt(obs.dates$start))[1]
      time.res <- difftime(as.POSIXlt(obs.dates$end[1]), as.POSIXlt(obs.dates$start[1]))
      start.hour <- as.POSIXlt(obs.dates$start)[1]$hour
      tz <- tryCatch(unlist(strsplit(obs.dates$start[1], split = "\\s"))[3], error = function(er) {er <- ""})
      if (is.null(names(sim.dates))) {
            tz.sim <- unlist(strsplit(sim.dates[[1]]$start[1], split = "\\s"))[3]
            sim.dates.ref <- as.POSIXlt(sim.dates[[1]]$start, tz = tz.sim)
      } else {
            tz.sim <- unlist(strsplit(sim.dates$start[1], split = "\\s"))[3]
            sim.dates.ref <- as.POSIXlt(sim.dates$start, tz = tz.sim)
      }
      aux.string <- paste(sim.dates.ref[1]$year + 1900, sim.dates.ref[1]$mon + 1, sim.dates.ref[1]$mday, start.hour, sep = "-")
      start <- seq(strptime(aux.string, "%Y-%m-%d-%H", tz), by = time.step, length.out = length(sim.dates.ref))
      end <- start + time.res
      usetz <- ifelse(identical(tz, ""), FALSE, TRUE)
      start <- format.POSIXct(start, "%Y-%m-%d %H:%M:%S", usetz = usetz)
      end <- format.POSIXct(end, "%Y-%m-%d %H:%M:%S", usetz = usetz)
      return(list("start" = start, "end" = end))
}
# End

