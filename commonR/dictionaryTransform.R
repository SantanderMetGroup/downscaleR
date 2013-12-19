# Description: Performs time bound definitions and variable transformation according to dictionary specifications
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dictionaryTransform <- function(Data, shortName, dictionary, timePars) {
	Data <- Data
	shortName <- shortName
	dictionary <- read.csv(dictionary)
	timePars <- timePars
	dicRow <- grep(shortName, dictionary$short_name)
	# Time bound definition
	ltb <- as.difftime(dictionary$lower_time_bound[dicRow], format = "%H", units = "hours")
	utb <- as.difftime(dictionary$upper_time_bound[dicRow], format = "%H", units = "hours")
	if (utb - ltb != 0) {
		dateSliceStart <- timePars$dateSlice - ltb
		dateSliceEnd <- timePars$dateSlice + utb
	} else {
		dateSliceStart <- timePars$dateSlice
		dateSliceEnd <- dateSliceStart
	}
	# Data transformation
	offset <- dictionary$offset[dicRow]
	scaleFac <- dictionary$scale[dicRow]
	Data <- Data * scaleFac + offset
	return(list("dateSliceStart" = dateSliceStart, "dateSliceEnd" = dateSliceEnd, "Data" = Data))
}
# End