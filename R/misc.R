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
      dimNames <- attr(obj$Data, "dimensions")
      aux <- if (any(grepl("var", dimNames))) {
            as.POSIXlt(obj$Dates[[1]]$start)$mon + 1      
      } else {
            as.POSIXlt(obj$Dates$start)$mon + 1      
      }
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
      dimNames <- attr(obj$Data, "dimensions")
      if (any(grepl("var", dimNames))) {
            aux.dates <- as.POSIXlt(obj$Dates[[1]]$start)
      } else {
            aux.dates <- as.POSIXlt(obj$Dates$start)
      }
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


#' @title Get geographical coordinates of a climate data object
#' @description Returns the coordinates of a climate data object, either stations
#'  or field
#' @param obj Any object extending the station or field classes
#' @return A list with x and y components
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @export
#' 

getCoordinates <- function(obj) {
      if ("station" %in% attr(obj$Data, "dimensions")) {
            x <- obj$xyCoords[ ,1]
            y <- obj$xyCoords[ ,2]
      } else {
            x <- obj$xyCoords$x
            y <- obj$xyCoords$y
      }
      return(list("x" = x, "y" = y))
}
# End



#' @title Set the 'dimensions' attribute 
#' @description Sets the 'dimensions' attribute of model out Data objects after downscaling
#' @param obs A observations object
#' @param multi.member Logical indicating if simulation data is a multimember
#' @return A character vector indicating the dimensions of the output object
#' @keywords internal
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @export

renameDims <- function(obs, multi.member) {
      dimNames <- attr(obs$Data, "dimensions")
      # Remove "station" from dimensions for single-station objects
      st.dim.index <- grep("station", dimNames)
      if (!identical(st.dim.index, integer(0))) {
            dim.st <- dim(obs$Data)[st.dim.index]
            if (identical(dim.st, 1L)) {
                  dimNames <- dimNames[-st.dim.index]
            }
      }
      if (isTRUE(multi.member)) {dimNames <- c("member", dimNames)}
      return(dimNames)
}
# End


#' Calculate the number of days of the current month
#' @param d A date (character) in format YYYY-MM-DD...
#' @return The number of days of the current month
#' @references 
#' \url{http://stackoverflow.com/questions/6243088/find-out-the-number-of-days-of-a-month-in-r}
#' @export

ndays <- function(d) {
      as.difftime(tail((28:31)[which(!is.na(as.Date(paste0(substr(d, 1, 8), 28:31), '%Y-%m-%d')))], 1), units = "days")
}
#End

#' Adjust time/start dates of a loaded object
#' @param timePars Object containing the relevant time parameters
#' @return A list with dates (POSIXct) start and end, defining the interval [start, end)
#' @details Sub-daily information is displayed only in case of subdaily data
#' @author J Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal

# timePars <- cube$timePars
adjustDates <- function(timePars) {
      interval <- 0
      if (timePars$aggr.m != "none") {
            mon.len <- sapply(timePars$dateSliceList, ndays)
            interval <- mon.len * 86400
      } else if (timePars$aggr.d != "none") {
            timePars$dateSliceList <- format(as.Date(substr(timePars$dateSliceList, 1, 10)), format = "%Y-%m-%d %H:%M:%S", usetz = TRUE) 
            interval <- 86400
      }
      formato <- ifelse(interval[1] == 0, "%Y-%m-%d %H:%M:%S", "%Y-%m-%d")
      dates.end <- format(as.POSIXct(as.POSIXlt(timePars$dateSliceList, tz = "GMT") + interval), format = formato, usetz = TRUE)
      dates.start <- format(as.POSIXct(as.POSIXlt(timePars$dateSliceList, tz = "GMT"), tz = "GMT"), format = formato, usetz = TRUE)
      return(list("start" = dates.start, "end" = dates.end))
}
# End

#' @title getIntersect
#' @description Get the common period of the objects obs and prd
#' @author S. Herrera
#' @export
#' @keywords internal

getIntersect <- function(obs,prd){
  dimNames <- attr(obs$Data, "dimensions")
  indDates <- which(as.POSIXct(obs$Dates$start, tz="GMT", format="%Y-%m-%d")==as.POSIXct(prd$Dates$start, tz="GMT", format="%Y-%m-%d"))
  auxDates <- as.POSIXct(obs$Dates$start[indDates], tz="GMT", format="%Y-%m-%d")
  indObs <- which(is.element(as.POSIXct(obs$Dates$start, tz="GMT", format="%Y-%m-%d"), auxDates))
  obs <- subsetDimension(obs, dimension = "time", indices = indObs)
  dimNames <- attr(prd$Data, "dimensions")
  indObs <- which(is.element(as.POSIXct(prd$Dates$start, tz="GMT", format="%Y-%m-%d"), auxDates))
  prd <- subsetDimension(prd, dimension = "time", indices = indObs)
  obj <- list(obs = obs, prd = prd)
  obj$Dates$start <- as.POSIXct(obs$Dates$start, tz="GMT", format="%Y-%m-%d")
  obj$Dates$end <- as.POSIXct(obs$Dates$end, tz="GMT", format="%Y-%m-%d")
  attr(obj$obs$Data, "dimensions") <- attr(obs$Data, "dimensions")
  attr(obj$prd$Data, "dimensions") <- attr(prd$Data, "dimensions")
  return(obj)
}



#' Identification of leap years
#' Identification of leap years
#' @param years a integer vector of (gregorian) years
#' @return a vector of indices of the position of leap years
#' @references \url{https://en.wikipedia.org/wiki/Leap_year}
#' @keywords internal
#' @export
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @examples
#' leap.years <- which.leap(1885:1937)
#' (1885:1937)[leap.years]

which.leap <- function(years) {
      which((years %% 4 == 0) & ((years %% 100 != 0) | years %% 400 == 0))
}


