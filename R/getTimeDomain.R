#' Selection of time slices of gridded datasets
#' 
#' Performs the selection of time slices based on season and year specification
#'  (including year-crossing seasons). Sub-routine of \code{loadGridDataset}
#'  
#' @param grid A java \sQuote{GeoGrid}
#' @param season Vector of months defining the requested season. Passed by 
#'  \code{loadGridDataset}.
#' @param years Vector of selected years. Passed by \code{loadGridDataset}.
#' @param verifTime Verification time defined as a numeric value (e.g. 6 for data
#' verifying at 06:00:00). Only applies for sub-daily datasets.
#' @returns A list of length two with the selected verification dates and a list
#'  with the index values defined as java objects of the class \sQuote{ucar.ma2.Range}.
#'  Output passed to \code{loadGridDataset}.
#' @details The indices of time positions are returned as a list, in order to read 
#' discontinuous data in the time dimension (i.e., seasons for different years)
#' from NetCDF files in a more efficient way. Note that this is not necessary 
#' in the case of ASCII source files, but implemented in all cases for simplicity
#' (in this case, time indices are unlisted)
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}

getTimeDomain <- function(grid, season, years, verifTime) {
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
    if (!identical(season, sort(season))) {
		if (years[1] == startYear) {
			warning(paste("First date in dataset: ", startDay, ". Seasonal data for the first year requested not available", sep = ""))
		} else {
			years <- append(years[1] - 1, years)
		}
		timeInd <- which((timeDates$year + 1900) %in% years & (timeDates$mon + 1) %in% season)
		ranks <- rank(season)
		for (i in 2:length(season)) {
			ranks[i] <- sign(ranks[i] - ranks[i-1])
		}
		crossSeason <- which(sign(ranks) == -1)
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
    if (!is.null(verifTime)) { 
        verifTimeInd <- which(dateSliceList[[1]]$hour == verifTime)
        if (length(verifTimeInd) == 0) {
            stop("Non-existing verification time selected.\nCheck value of argument 'verifTime'")
        }
        timeShift <- as.integer(-(verifTimeInd[1] - 1))
        for (i in 1:length(dateSliceList)) {
            dateSliceList[[i]] <- dateSliceList[[i]][verifTimeInd]
        }
    }
    if (is.null(verifTime)) {
        timeStride <- 1L
        timeShift <- 0L
    } else {
        timeStride <- as.integer(verifTimeInd[2] - verifTimeInd[1])
    }
    dateSlice <- do.call("c", dateSliceList)
    tRanges <- lapply(1:length(timeIndList), function(j) .jnew("ucar/ma2/Range", as.integer(timeIndList[[j]][1]), as.integer(timeIndList[[j]][length(timeIndList[[j]])]), timeStride)$shiftOrigin(timeShift))
    return(list("dateSlice" = dateSlice, "timeResInSeconds" = timeResInSeconds, "tRanges" = tRanges))
}
# End