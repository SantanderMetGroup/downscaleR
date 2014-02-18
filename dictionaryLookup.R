# Description: Finds the standard variable in the dictionary and returns the standardName
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dictionaryLookup <- function(dictionary, var) {
      dictionary <- tryCatch({read.csv(dictionary)}, error = function(e) stop("Dictionary not found"))
      dicRow <- grep(paste("^", var, "$", sep = ""), dictionary$identifier)
      if (length(dicRow) == 0) {
            stop("Variable requested does not match any identifier in the dictionary")
      } else {
            shortName <- as.character(dictionary$short_name[dicRow])
						standardName <- as.character(dictionary$identifier[dicRow])
      }
      return(list("ShortName" = shortName, "StandardName" = standardName))
}
# End
