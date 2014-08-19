#' Time  index positions for station dataset selections
#' 
#' Get time index positions for loading ascii station data
#' 
#' @param timeDates a POSIXlt vector of time dates
#' @param season A vector of months defining the season selected
#' @param years A vector of (continuous) year selection
#' @return A list with a vector of time index positions and the corresponding POSIXlt dates
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal

getTimeDomainStations <- function(timeDates, season, years) {
    if (is.null(season)) {
        season <- 1:12
    }
    allYears <- unique(timeDates$year + 1900)
    startYear <- head(allYears, 1L)
    endYear <- tail(allYears, 1L)
    if (is.null(years)) {
        years <- allYears
    } 
    if (years[1] < startYear & tail(years, 1L) > endYear) {
        warning("Year selection out of dataset range. Only available years will be returned")
        years <- allYears
    }
    if (years[1] < startYear) {
        warning("First year in dataset: ", startYear,". Only available years will be returned")
        years <- startYear:years[length(years)]
    }
    if (tail(years, 1L) > endYear) {
        warning("Last year in dataset: ", endYear,". Only available years will be returned")
        years <- years[1]:endYear
    }
    # Year-crossing seasons - year to take the initialization
    if (!identical(season, sort(season))) {
          if (years[1] == startYear) { 
                warning(paste("First forecast day in dataset: ", timeDates[1], ".\nRequested seasonal data for ", startYear," not available", sep=""))
                years <- years[-length(years)]
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
    }else{
          timeInd <- which((timeDates$year + 1900) %in% years & (timeDates$mon + 1) %in% season)
    }  
    timeDates <- timeDates[timeInd]
    return(list("timeInd" = timeInd, "timeDates" = timeDates))
}
