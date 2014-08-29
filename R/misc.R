#'@title Get season from a station or field object
#'@description Retrieves the season encompassed by a station or field object
#'@param obj Any object extending the station or field classes
#'@return An integer vector with the season
#'@keywords internal
#'@author J. Bedia
#'@export

getSeason <- function(obj) {
      aux <- as.POSIXlt(obj$Dates$start)$mon + 1      
      return(unique(aux))
}
# End