#' Selection of time slices of gridded datasets
#' 
#' Performs the selection of time slices based on season and year specification
#'  (including year-crossing seasons). Sub-routine of \code{loadGridDataset}
#'  
#' @param grid A java \sQuote{GeoGrid}
#' @param season Vector of months defining the requested season. Passed by 
#'  \code{loadGridDataset}.
#' @param years Vector of selected years. Passed by \code{loadGridDataset}.
#' @param time Verification time defined as a character string (e.g. \dQuote{06} for data
#' verifying at 06:00:00). Only applies for sub-daily datasets.
#' @returns A list of length two with the selected verification dates and a list
#'  with the index values defined as java objects of the class \sQuote{ucar.ma2.Range}.
#'  Output passed to \code{loadGridDataset}.
#' @details The indices of time positions are returned as a list, in order to read 
#' discontinuous data in the time dimension from NetCDF files in a more efficient way
#' (i.e. seasons for different years). 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}

getTimeDomain <- function(grid, dic, season, years, time) {
    message("[", Sys.time(), "] Defining time selection parameters")
    gcs <- grid$getCoordinateSystem()
    timeDates <- javaCalendarDate2rPOSIXlt(gcs$getTimeAxis()$getCalendarDates())
    timeResInSeconds <- gcs$getTimeAxis()$getTimeResolution()$getValueInSeconds()
    startDay <- timeDates[1]
    endDay <- timeDates[length(timeDates)]
	startYear <- startDay$year + 1900
	endYear <- endDay$year + 1900
    if (is.null(years)) {
	      years <- as.integer(startYear : endYear)
	}
	if (years[1] < startYear | years[length(years)] > endYear) {
		warning("Year selection out of boundaries. Only available years will be returned")
	}
	if (years[1] < startYear) {
		years <- startYear : years[length(years)]
	}
    if (years[length(years)] > endYear) {
	    years <- years[1] : endYear
	}
	if (is.null(season)) {
		season <- as.integer(1:12)
	} else {
		season <- as.integer(season)
		if (min(season) < 1 | max(season) > 12) {
			stop("Invalid season definition")
		}
	}
    if (!isTRUE(identical(season, sort(season)))) {
	    if (years[1] == startYear) {
			warning(paste("First date in dataset: ", startDay, ". Seasonal data for the first year requested not available", sep = ""))
		} else {
			years <- append(years[1] - 1, years)
		}
		timeInd <- which((timeDates$year + 1900) %in% years & (timeDates$mon + 1) %in% season)
		crossSeason <- which(c(1, diff(season)) < 0)
            rm.ind <- which((timeDates$mon + 1) %in% season[1 : (crossSeason - 1)] & (timeDates$year + 1900) %in% years[length(years)])
		if (length(years) > 1) {
			rm.ind <- c(rm.ind, which((timeDates$mon + 1) %in% season[crossSeason : length(season)] & (timeDates$year + 1900) %in% years[1]))
		}
		timeInd <- setdiff(timeInd, rm.ind)
    } else {
		timeInd <- which((timeDates$year + 1900) %in% years & (timeDates$mon + 1) %in% season)
	}
    dateSlice <- timeDates[timeInd]
    timeIndList <- list()
    dateSliceList <- list()
    if (length(dateSlice) > 1) {
        brkInd <- rep(1, length(timeInd))
            for (i in 2:length(timeInd)) {
		      brkInd[i] <- timeInd[i] - timeInd[i-1]
	      }
	      brkInd <- c(1, which(brkInd > 1), length(timeInd) + 1)
            if (length(brkInd) == 0) { 
		      timeIndList[[1]] <- timeInd - 1
                  dateSliceList[[1]] <- dateSlice
	        } else {
		    for (i in 2:length(brkInd)) {
		        timeIndList[[i - 1]] <- timeInd[brkInd[i - 1] : (brkInd[i] - 1)] - 1
                dateSliceList[[i - 1]] <- dateSlice[brkInd[i - 1] : (brkInd[i] - 1)]
		    }
	    }
	} else {
	    timeIndList[[1]] <- timeInd - 1
	    dateSliceList[[1]] <- dateSlice
    }
      # Sub-routine for setting stride and shift along time dimension    
      if (is.null(dic) & time != "none") {
            stop("Time resolution especification incompatible with non-standard variable requests\nUse the dictionary or set the 'time' argument to NULL")
      }
      if (is.null(dic) | isTRUE(dic$doDailyMean)) {
            timeStride <- 1L
            timeShift <- 0L
      } else {
            if (time == "DD" | time == "none") {
                  timeStride <- 1L
                  timeShift <- 0L
            } else {
                  time <- as.integer(time)
                  timeIndList <- lapply(1:length(dateSliceList), function(x) {
                        which(dateSliceList[[x]]$hour == time)
                  })
                  if (length(timeIndList[[1]]) == 0) {
                        stop("Non-existing verification time selected.\nCheck value of argument 'time'")
                  }
                  dateSliceList <- lapply(1:length(dateSliceList), function(x) {
                        dateSliceList[[x]][timeIndList[[x]]]
                  })
                  timeStride <- as.integer(diff(timeIndList[[1]])[1])
                  timeShift <- as.integer(-(timeIndList[[1]][1] - 1))
                  timeIndList <- NULL
            }
      }
      dateSlice <- do.call("c", dateSliceList)
      dateSliceList <- NULL
    # Sub-routine for adjusting times in case of deaccumulation
    deaccumFromFirst <- NULL
    if (!is.null(dic)) {
        if (dic$deaccum == 1) {
              if (timeIndList[[1]][1] > 1) {
                    deaccumFromFirst <- FALSE
                    timeIndList <- lapply(1:length(timeIndList), function(x) {
                          c(timeIndList[[x]][1] - 1, timeIndList[[x]])
                    })
              } else {
                    deaccumFromFirst <- TRUE
              }
        }
    }
    # Sub-routine for calculation of time bounds
    dateSlice <- timeBounds(dic, dateSlice)
    tRanges <- lapply(1:length(timeIndList), function(j) .jnew("ucar/ma2/Range", as.integer(timeIndList[[j]][1]), as.integer(tail(timeIndList[[j]], 1L)), timeStride)$shiftOrigin(timeShift))
    return(list("dateSlice" = dateSlice, "timeResInSeconds" = timeResInSeconds, "tRanges" = tRanges))
}
# End