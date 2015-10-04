#' @title Select an arbitrary subset from a field or multifield along one of its dimensions
#'
#' @description Create a new field/multifield that is a subset of the input field 
#' along the selected dimension
#'
#' @param field The input field to be subset. This is either a field, as returned e.g. by \code{loadGridData}, a
#' multifield, as returned by \code{loadMultiField} or \code{makeMultiField}, or other types of multimember fields
#' (possibly multimember multifields) as returned e.g. by \code{\link[ecomsUDG.Raccess]{loadECOMS}}.
#' @param dimension Character vector indicating the dimension along which the positions indicated by the \code{indices} paraneter.
#' @param indices An integer vector indicating \strong{the positions} of the dimension to be extracted.
#' @return A new field object that is a logical subset of the input field along the specified dimension.
#' @details
#' 
#' The attribute \code{subset} will be added taking the value of the \code{dimension} parameter.
#'  
#' @importFrom abind asub
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com} and S. Herrera
#' @export
#' @family subsetting
#' @examples
#' # Example - Member subset
#' data(tasmax_forecast)
#' plotMeanField(tasmax_forecast, TRUE)
#' # Selection of a smaller domain over the Iberian Peninsula and members 3 and 7
#' sub <- subsetDimension(tasmax_forecast,
#'                    dimension = "member",
#'                    indices = c(1,3))
#' plotMeanField(sub, multi.member = TRUE)

subsetDimension <- function(field, dimension = NULL, indices = NULL) {
  dimNames <- attr(field$Data, "dimensions")
  if (!is.null(indices)){
    field$Data <- asub(field$Data, indices, grep(dimension, dimNames))
    attr(field$Data, "dimensions") <- dimNames
    if ("time" %in% dimension){
      field$Dates$start <- field$Dates$start[indices]
      field$Dates$end <- field$Dates$end[indices]
    }
    if ("lon" %in% dimension){
      field$xyCoords$x <- field$xyCoords$x[indices]
    }
    if ("lat" %in% dimension){
      field$xyCoords$y <- field$xyCoords$y[indices]
    }
    if ("member" %in% dimension){
      field$Members <- field$Members[indices]
      if (is.list(field$InitializationDates)) { # e.g. CFSv2 (members defined through lagged runtimes)
        field$InitializationDates <- field$InitializationDates[indices]
      } 
    }
    attr(field$Variable, "subset") <- dimension
  }else{
    warning("Argument 'indices' is NULL and none subsetting has been applied. The same input 'field' is returned.")
  }
  return(field)
}
# End
