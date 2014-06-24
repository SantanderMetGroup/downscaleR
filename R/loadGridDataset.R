#' Loads a user-defined subset of a gridded CDM dataset
#' 
#' Loads a user-defined subset from a gridded dataset compliant with the Common
#'  Data Model interface
#' 
#' @param dataset A complete URL to the file describing the dataset (usually a NcML)
#' @param var A character string indicating the variable requested. This is usually
#' a variable defined in the vocabulary. See details.
#' @param vocabulary Type of vocabulary to be used for data homogeneization. Currently
#' only two options are accepted: \dQuote{\code{standard}} (the default) or \dQuote{\code{none}}.
#'  The latter is used only for original (non-homogenized) data retrieval. See details.
#' @param dictionary Full path to the dictionary file (.dic). Default to \code{NULL},
#'  indicating that the .dic file is in the same directory as the \code{dataset}. See details.
#' @param lonLim
#' @param latLim
#' @param level A numeric value indicating the vertical level.
#' @param season
#' @param years 
#' @return
#' @export
#' @details
#' @seealso \code{\link{dataInventory}}
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @references \url{http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/tutorial/GridDatatype.html}

loadGridDataset <- function(var, grid, dic, level, season, years, time, latLon) {
    timePars <- getTimeDomain(grid, dic, season, years, time)
    levelPars <- getVerticalLevelPars(grid, level)
    mdArray <- makeSubset(grid, timePars$tRanges, levelPars$zRange, latLon)
    if (!is.null(dic)) {
        isStandard <- TRUE
        mdArray <- dictionaryTransformGrid(dic, timePars, mdArray)
    } else {
        isStandard <- FALSE
    }
    if (isTRUE(latLon$revLat)) {
        mdArray <- revArrayLatDim(mdArray, grid)
    }
    return(list("Variable" = list("varName" = var, "isStandard" = isStandard),
                "Data" = mdArray,
                "xyCoords" = c(latLon$xyCoords, "CRS_string" = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"), 
                "Dates" = timePars$dateSlice))
}
# End

 


