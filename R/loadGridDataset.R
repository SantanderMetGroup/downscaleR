#' Loads a user-defined subset of a gridded CDM dataset
#' 
#' Loads a user-defined subset from a gridded dataset compliant with the Common
#'  Data Model interface
#' 
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @export
#' @keywords internal

loadGridDataset <- function(var, grid, dic, level, season, years, time, latLon, aggr.d, aggr.m) {
      timePars <- getTimeDomain(grid, dic, season, years, time, aggr.d, aggr.m)
      levelPars <- getVerticalLevelPars(grid, level)
      cube <- makeSubset(grid, timePars, levelPars, latLon)
      timePars <- NULL
      if (!is.null(dic)) {
            isStandard <- TRUE
            cube$mdArray <- dictionaryTransformGrid(dic, cube$timePars, cube$mdArray)
      } else {
            isStandard <- FALSE
      }
      if (isTRUE(latLon$revLat)) {
            cube$mdArray <- revArrayLatDim(cube$mdArray, grid)
      }
      Variable <- list("varName" = var, "level" = levelPars$level)
      attr(Variable, "is_standard") <- isStandard
      if (isStandard) {
            data(vocabulary, envir = environment())
            attr(Variable, "units") <- as.character(vocabulary[grep(paste0("^", var, "$"), vocabulary$identifier,), 3])
            attr(Variable, "longname") <- as.character(vocabulary[grep(paste0("^", var, "$"), vocabulary$identifier,), 2])
      } else {
            attr(Variable, "units") <- "undefined"
            attr(Variable, "longname") <- "undefined"
      }
      attr(Variable, "daily_agg_cellfun") <- cube$timePars$aggr.d
      attr(Variable, "monthly_agg_cellfun") <- cube$timePars$aggr.m
      attr(Variable, "verification_time") <- time
      out <- list("Variable" = Variable, "Data" = cube$mdArray, "xyCoords" = latLon$xyCoords, "Dates" = adjustDates(cube$timePars))
      return(out)
      
}
# End
    
    


