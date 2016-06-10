#' @title Get season 
#' @description Retrieves the season encompassed by a station or grid object
#' @param obj Any object extending the station or grid classes
#' @return An integer vector with the season
#' @author J. Bedia 
#' @export
#' @examples 
#' data(iberia_ncep_ta850)
#' getSeason(iberia_ncep_ta850) # Boreal winter (DJF)

getSeason <- function(obj) {
      if ("season" %in% names(attributes(obj$Dates))) {
            attr(obj$Dates, "season")
      } else {
            dimNames <- attr(obj$Data, "dimensions")
            aux <- if (any(grepl("var", dimNames))) {
                  as.POSIXlt(obj$Dates[[1]]$start)$mon + 1      
            } else {
                  as.POSIXlt(obj$Dates$start)$mon + 1      
            }
            unique(aux)
      }
}
# End

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

#' @title Get years as a factor
#' @description Extract the year as a factor (e.g. for computing annual statistics)
#' @param obj Any object extending the station or grid classes
#' @return A vector of years of the same length as the time dimension of the object, 
#' seasonally-adjusted in the case of year-crossing seasons (e.g. DJF). See details.
#' @details The function performs a very basic operation, extracting the year element from the 
#' dates previously converted to POSIXlt. The trick lies in the year-crossing seasons. For instance:
#'  by convention, winter 2001 encompasses December 2000 and January, February 2001. Therefore, in order to compute
#' annual statistics for a year-crossing season, it is necessary to modify first the vector of years, 
#' and assign year 2001 to the preceding December. Similarly, the next December 2001 belongs to winter 2002,
#'  and so on... The function is useful for computing and/or plotting annual statistics, seasonal climatologies ... 
#' @note Warning:
#' The function should no be used to extract the vector of actual date years
#' @author J. Bedia 
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
      aux.dates <- if (any(grepl("var", dimNames))) {
            obj$Dates[[1]]$start
      } else {
            obj$Dates$start
      }
      yrs <- as.numeric(substr(aux.dates,1,4))
      mon <- as.numeric(substr(aux.dates,6,7))
      if (identical(yrs, unique(yrs))) {
          yrs
      } else {    
          if (!identical(season, sort(season))) {
              yy <- unique(yrs)[-1]
              aux <- match(mon, season)
              brks <- c(1, which(diff(aux) < 0) + 1, length(aux) + 1)
              l <- lapply(1:(length(brks) - 1), function(x) {
                  a <- yrs[brks[x]:(brks[x + 1] - 1)]
                  return(rep(yy[x], length(a)))
              })
              yrs  <- do.call("c", l)
          }
      }
      return(yrs)
}
# End



#' @title Set the 'dimensions' attribute 
#' @description Sets the 'dimensions' attribute of model out Data objects after downscaling
#' @param obs A observations object
#' @param multi.member Logical indicating if simulation data is a multimember
#' @return A character vector indicating the dimensions of the output object
#' @keywords internal
#' @author J. Bedia


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


#' @title getIntersect
#' @description Get the common period of the objects obs and prd
#' @author S. Herrera
#' @keywords internal

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



#' Identification of leap years
#' Identification of leap years
#' @param years a integer vector of (gregorian) years
#' @return a vector of indices of the position of leap years
#' @references \url{https://en.wikipedia.org/wiki/Leap_year}
#' @keywords internal
#' @export
#' @author J. Bedia 
#' @examples
#' leap.years <- which.leap(1885:1937)
#' (1885:1937)[leap.years]

which.leap <- function(years) {
      which((years %% 4 == 0) & ((years %% 100 != 0) | years %% 400 == 0))
}


#' @title Land borders
#' @description Add land borders to a map
#' @param ... Graphical parameters passed to \code{\link{lines}}.
#' @return Draws a simplied land border areas as lines onto the map
#' @details The function loads a built-in world segments dataset created ad hoc to avoid dependencies on other packages (i.e. 'maps').
#' Geographical lonlat coordinates in wgs84.
#' @source Postprocessed from the original shapefile from Natural Earth (http://www.naturalearthdata.com/downloads/110m-physical-vectors/)
#' @author J. Bedia
#' @keywords internal
#' @export

draw.world.lines <- function(...) {
      load(file.path(find.package("downscaleR"), "wrl.Rda"))
      for (i in 1:length(node.list)) {
            lines(node.list[[i]][,1], node.list[[i]][,2], ...)            
      }
}

#' @title  Retrieve dimensions attribute
#' @description Retrieve dimensions attribute
#' @param obj A grid or station object
#' @return A character vector with the dimensions attribute of the object's \code{Data} component.
#' @keywords internal
#' @author J. Bedia

getDim <- function(obj) {
      attr(obj[["Data"]], "dimensions")
}

#' @title  Retrieve array shape 
#' @description Retrieve array attributes 'dimensions' and 'dim'
#' @param obj A grid or station object
#' @return An integer vector with dim values, labelled with the \code{"dimension"} attribute names
#' @keywords internal
#' @author J. Bedia

getShape <- function(obj, dimension = NULL) {
      dimNames <- getDim(obj)
      shape <- dim(obj[["Data"]])
      if (!is.null(dimension)) {
            ind <- match(dimension, dimNames)
            if (anyNA(ind)) stop("Input 'dimension' value not found")
            shape <- shape[ind]
            dimNames <- dimNames[ind]
      }
      names(shape) <- dimNames
      return(shape)
}


