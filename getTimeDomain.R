# Author: juaco
# Date: 28 Oct 2013
# Info: netcdf-Java 4.3.18 / R 3.0.2/ rJava 0.9-4
# TODO: complete documentation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ====================
# getTimeDomain.R #
# ====================
# Description: Function to perform the selection of time slices based on season and year specification (including year-crossing seasons)
# Arguments: # 'timeDates': a POSIXlt vector containing all dates in the dataset
		     # 'season': arguments passed by calling function
		     # 'years': arguments passed by calling function
# Value: A list with two elements:
		   # [[1]] 'dateSlice': A POSIXlt vector of selected dates
		   # [[2]] 'timeIndList' A list with index positions of selected dates. Each element of the list corresponds to a continuous chunk of dates. See details.
# Details: The indices of time positions are returned as a list, in order to read discontinuous data in the time dimension
			# (i.e., seasons for different years) from NetCDF files in a more efficient way.
			# Note that this is not necessary in the case of ASCII source files, but implemented in all cases for simplicity (in this case, time indices are unlisted)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

getTimeDomain <- function(timeDates, season, years) {
	timeDates <- timeDates
	startDay <- timeDates[1]
	endDay <- timeDates[length(timeDates)]
	startYear <- startDay$year + 1900
	endYear <- endDay$year + 1900
	if (is.null(years)) {
		years <- startYear : endYear
	}
	years <- as.integer(years)
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
	## year-crossing season
	if (identical(season, sort(season)) == FALSE) {
		if (years[1] == startYear) {
			warning(paste("First forecast date in dataset: ", startDay, ". Seasonal data for the first year requested not available", sep = ""))
		} else {
			years <- append(years[1] - 1, years)
		}
		timeInd <- which((timeDates$year + 1900) %in% years & (timeDates$mon + 1) %in% season)
		# Remove tails/heads of the years to match season start/end
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
	# Particion de la time slice en tramos continuos en forma de lista (timeIndList)
	brkInd <- rep(1, length(timeInd))
	for (i in 2:length(timeInd)) {
		brkInd[i] <- timeInd[i] - timeInd[i-1]
	}
	brkInd <- c(1, which(brkInd > 1), length(timeInd) + 1)
	timeIndList <- list()
	if (length(brkInd) == 0) { # Si no hay discontinuidades
		timeIndList[[1]] <- timeInd - 1
	} else {
		for (i in 2:length(brkInd)) {
			timeIndList[[i - 1]] <- timeInd[brkInd[i - 1] : (brkInd[i] - 1)] - 1
		}
	}
	return(list("dateSlice" = dateSlice, "timeIndList" = timeIndList))
}
# End
