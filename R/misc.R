#' @title Get season from a station or field object
#' @description Retrieves the season encompassed by a station or field object
#' @param obj Any object extending the station or field classes
#' @return An integer vector with the season
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @export
#' @examples 
#' data(iberia_ncep_ta850)
#' getSeason(iberia_ncep_ta850) # Boreal winter (DJF)
#' 

getSeason <- function(obj) {
      aux <- as.POSIXlt(obj$Dates$start)$mon + 1      
      return(unique(aux))
}
# End

#' @title Compute dates of a downscaled observational dataset
#' @description The function calculates the appropiate Dates slot in the returned output of downscaling functions,
#' considering the possible mismatches in time resolution between predictors and predictand, the multifield dates slot etc.
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

#' @title Get years as a factor
#' @description Extract the year as a factor (e.g. for computing annual statistics)
#' @param obj Any object extending the station or field classes
#' @return A vector of years of the same length as the time dimension of the object, 
#' seasonally-adjusted in the case of year-crossing seasons (e.g. DJF). See details.
#' @details The function performs a very basic operation, extracting the year element from the 
#' dates previously converted to POSIXlt. The trick lies in the year-crossing seasons. For instance:
#'  by convention, winter 2001 encompasses December 2000 and January, February 2001. Therefore, in order to compute
#' annual statistics for a year-crossing season, it is necessary to modify first the vector of years, 
#' and assign year 2001 to the preceding December. Similarly, the next December 2001 belongs to winter 2002,
#'  and so on... The function is useful for computing and/or plotting annual statistics, seasonal climatologies ... 
#' @section Warning:
#' The function should no be used to extract the actual years vector
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @export
#' @examples 
#' data(iberia_ncep_hus850)
#' getSeason(iberia_ncep_hus850)
#' # Winter 1991-2010
#' range(iberia_ncep_hus850$Dates$start)
#' ## Time series for the first point
#' # Dates vector
#' time <- as.POSIXlt(iberia_ncep_hus850$Dates$start, tz = "GMT")
#' hus850 <- iberia_ncep_hus850$Data[ ,1,1]
#' plot(time, hus850, ty = "l")
#' ## Computation of the annual series for winter specific humidity:
#' par(mfrow = c(2,1))
#' ## Wrong:
#' years <- as.POSIXlt(iberia_ncep_hus850$Dates$start)$year + 1900
#' x <- tapply(hus850, INDEX = list(years), FUN = mean)
#' plot(unique(years), x, ty = "b")
#' points(1990, x[1], col = "red", cex = 2, lwd = 2)
#' ## Correct:
#' years <- getYearsAsINDEX(iberia_ncep_hus850)
#' x <- tapply(hus850, INDEX = years, FUN = mean)
#' plot(unique(years), x, ty = "b")
#' par(mfrow = c(1,1))
#' 

getYearsAsINDEX <- function(obj) {
      season <- getSeason(obj)
      aux.dates <- as.POSIXlt(obj$Dates$start)
      yrs <- aux.dates$year + 1900
      if (!identical(season, sort(season))) {
            yy <- unique(yrs)[-1]
            aux <- match(aux.dates$mon + 1, season)
            brks <- c(1, which(diff(aux) < 0) + 1, length(aux) + 1)
            l <- lapply(1:(length(brks) - 1), function(x) {
                  a <- yrs[brks[x] : (brks[x + 1] - 1)]
                  return(rep(yy[x], length(a)))
            })
            yrs  <- do.call("c", l)
      }
      return(yrs)
}
# End
