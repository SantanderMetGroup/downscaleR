#' @title Searches variable string in the dictionary
#' 
#' @description Searches variable string provided in the dictionary to map it in the vocabulary, in order to get
#' all the necessary information for variable homogenization. It also includes a new column specifying the
#' aggregation function to be applied (if any).
#' 
#' @param dicPath Full path to the dictionary file (a csv file with extension \sQuote{.dic}).
#' @param var Character string with the (standard) name of the variable
#' @param time Time specification.
#' 
#' @return A data.frame of 1 row with the mapping information
#' 
#' @references \url{http://meteo.unican.es/ecoms-udg/RPackage/Homogeneization}
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' 
#' @keywords internal

dictionaryLookup <- function(dicPath, var, time) {
      message("[", Sys.time(), "] Defining homogeneization parameters for variable \"", var, "\"")
      dictionary <- tryCatch({read.csv(dicPath, stringsAsFactors = FALSE)}, error = function(e) stop("Dictionary not found"))
      dicRow <- grep(paste("^", var, "$", sep = ""), dictionary$identifier) 
      if (length(dicRow) == 0) {
            stop("Variable requested does not match any identifier in the dictionary\nType 'help(vocabulary)' for help on standard variable naming")
      }
      dailyAggr <- NA
      if (length(dicRow) > 1) {
            if (time == "DD") {
                  dicRow <- dicRow[dictionary$time_step[dicRow] == "24h"]
                  if (length(dicRow) == 0) {
                        dicRow <- grep(paste("^", var, "$", sep = ""), dictionary$identifier)                  
                        dicRow <- dicRow[dictionary$time_step[dicRow] == "6h"]
                  }
            } else {
                  dicRow <- dicRow[dictionary$time_step[dicRow] == "6h"]
            }
      } else {
            if (dictionary$time_step[dicRow] == "12h" & time == "DD") {
                  stop("Cannot compute daily mean from 12-h data")
            }
            if ((time == "06" | time == "18") & dictionary$time_step[dicRow] == "12h") {
                  stop("Requested 'time' value (\"", time, "\") not available for 12-h data")
            }
            if ((time != "none" & time != "DD") & (dictionary$time_step[dicRow] == "24h")) {
                  stop("Subdaily data not available for variable \"", var, "\". Check value of argument 'time'")
            }
            if (time == "DD" & dictionary$time_step[dicRow] == "24h") {
                  time <- "none"
            }
            if (time == "DD") {
                  dailyAggr <- "mean"
                  if (var == "tp" | var == "rlds" | var == "rsds") {
                        dailyAggr <- "sum"
                        message("NOTE: daily accumulated will be calculated from the 6-h model output")
                  } else {
                        message("NOTE: daily mean will be calculated from the 6-h model output")
                  }
            }
      }
      dic <- cbind.data.frame(dictionary[dicRow, ], "dailyAggr" = I(dailyAggr))
      return(dic)
}
# End

