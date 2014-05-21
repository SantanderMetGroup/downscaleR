#' Performs conversion between non-standard (dataset) and standard (vocabulary) variables
#' 
#' Finds the standard variable in the vocabulary matching the corresponding row
#'  in the dictionary, containing all the necessary information for variable
#'  transformation. Sub-routine of \code{loadGridDataset}
#'  
#'  @param dictionary A valid URL to the .dic file. Argument passed by the loading function
#'  @param var Character string indicating the variable requested. Argument passed by the loading function
#'  @return A data frame with one row, whose columns correspond to the different fields of
#'  the dictionary
#'  @author J. Bedia \email{joaquin.bedia@@gmail.com}

dictionaryLookup <- function(dictionary, var) {
    dictionary <- tryCatch({read.csv(dictionary, stringsAsFactors = FALSE)}, error = function(e) stop("Dictionary not found"))
    dicRow <- grep(paste("^", var, "$", sep = ""), dictionary$identifier)
    if (length(dicRow) == 0) {
        stop("Variable requested does not match any identifier in the dictionary")
    }
    return(as.data.frame(dictionary[dicRow, ]))
}
# End



