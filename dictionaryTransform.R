# Description: Performs time bound definitions and variable transformation according to dictionary specifications
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dictionaryTransform <- function(Data, dic, timePars) {
	ltb <- as.difftime(dic$lower_time_bound, format = "%H", units = "hours")
	utb <- as.difftime(dic$upper_time_bound, format = "%H", units = "hours")
	
      if (utb - ltb != 0) {
		dateSliceStart <- timePars$dateSlice - ltb
		dateSliceEnd <- timePars$dateSlice + utb
	} else {
		dateSliceStart <- timePars$dateSlice
		dateSliceEnd <- dateSliceStart
	}
	Data <- Data * dic$scale + dic$offset
	return(list("dateSliceStart" = dateSliceStart, "dateSliceEnd" = dateSliceEnd, "Data" = Data))
}
# End