# Author: juaco
# Date: 28 Oct 2013
# Info: netcdf-Java 4.3.18 / R 3.0.2/ rJava 0.9-4
# TODO: complete documentation
# TODO: implement netCDF version
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ====================
# loadObservations.R #
# ====================
# Description: loads observations data from station datasets.
# Arguments: # 'source.dir': A valid path to the directory containing the station files
			 # 'file.format': Wether the stations data are stored in a netCDF or ASCII (default) file. See details for standard format definition.
			 # 'var': character string indicating the variable to be loaded. Note that the notation depends on wether the dictionary is being used or not.
			 # 'dictionary': Optional. A full path to the file containing the dictionary. See details
			 # 'stationID': Optional. A character vector indicating the code names of the stations to be loaded.
			 # 'lonLim':
			 # 'latLim':
			 # 'season':
			 # 'years':
# Value: a list with the following elements:
             # 'VarName': name of the variable requested (Note that the name can be either standard or not depending on wether a dictionary is being used)
			 # 'isStandard': logical flag indicating if the variable is standard (dictionary is in use) or not.
			 # 'Stations': A data frame containing the station codes and their names (locations).
             # 'LonLatCoords': A 2-D matrix with longitude and latitudes of the stations
			 # 'Dates': A list with the verification times of each record in the time series. This is represented by a list with two elements:
					# 'Start': The lower time bound of the verification period
					# 'End': The upper time bound of the verification period
			 # 'Data': A matrix with the values of the variable. Dates are ordered by rows and Stations in columns.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
loadObservations <- function(source.dir, file.format = c("ascii", "netcdf"), var, dictionary = NULL, stationID = NULL, lonLim = NULL, latLim = NULL, season = NULL, years = NULL) {
	source.dir <- source.dir
    ff <- match.arg(file.format, c("ascii", "netcdf"))
    var <- var
	dictionary <- dictionary
    stationID <- stationID 
    lonLim <- lonLim
    latLim <- latLim
    season <- season
    years <- years
	if (is.null(dictionary)) {
		shortName <- var
		standardName <- var
		isStandard <- FALSE
	} else {
		shortName <- dictionaryLookup(dictionary, var)$ShortName
		standardName <- dictionaryLookup(dictionary, var)$StandardName
		isStandard <- TRUE
	}
    if (ff == "ascii") {
    	out <- loadObservations.ASCII(source.dir, shortName, stationID, lonLim, latLim, season, years)
    } else if (ff == "netcdf") {
        out <- loadObservations.NetCDF(source.dir, shortName, stationID, lonLim, latLim, season, years)
    }
	if (is.null(dictionary) == FALSE) {
		trans <- dictionaryTransform(out$Data, shortName, dictionary, out$TimePars)
		out$Data <- trans$Data
		dateSliceStart <- trans$dateSliceStart
		dateSliceEnd <- trans$dateSliceEnd
	} else {
		varTimeStep <- difftime(out$TimePars$dateSlice[2], out$TimePars$dateSlice[1])
	  	dateSliceStart <- out$TimePars$dateSlice
	  	dateSliceEnd <- out$TimePars$dateSlice + varTimeStep
  	}
	dateSliceStart <- as.POSIXlt(dateSliceStart)
	dateSliceEnd <- as.POSIXlt(dateSliceEnd)
	return(list("VarName" = standardName, "isStandard" = isStandard, "Stations" = out$Stations, "LonLatCoords" = out$LonLatCoords,  "Dates" = list("Start" = dateSliceStart, "End" = dateSliceEnd), "Data" = out$Data))
}
# End      
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Example
# list.files('/home/juaco/temp/VALUE_validationTest_v2')
# source.dir = '/home/juaco/temp/VALUE_validationTest_v2'
#dictionary = "/home/juaco/temp/VALUE_validationTest_v2/dictionary.dic"
#var = "pr"
# dictionary = NULL
# var = "temp"
# lonLim = c(-10,10)
# latLim = c(35,44)
# ex <- loadObservations(source.dir, var = var, stationID = NULL, lonLim = c(-10,10), latLim = c(35,44), dictionary = dictionary, season = c(12,1,2), years = 1987:2000)
#str(ex)
#ex <- loadObservations(source.dir, var = var, stationID = c('000420','001394'), dictionary = dictionary, season = c(12,1,2), years = 1987:2000)	  







