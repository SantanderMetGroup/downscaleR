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

loadGridDataset <- function(dataset, var, vocabulary = "standard", dictionary = NULL, lonLim = NULL, latLim = NULL, level = NULL, season = NULL, years = NULL, verifTime = NULL) {
    vocabulary <- match.arg(vocabulary, choices = c("standard", "none"))
    if (vocabulary == "standard") {
        if (is.null(dictionary)) {
            dictionary <- gsub("ncml$|nc$", "dic", dataset)
        } else {
            dictionary <- dictionary
        }
        dic <- dictionaryLookup(dictionary, var)
        shortName <- dic$short_name
        isStandard <- TRUE
    } else {
        warning("Non-standard variable requested: \'", var, "\'")
        shortName <- var
        isStandard <- FALSE
    }
    gds <- J("ucar.nc2.dt.grid.GridDataset")$open(dataset)
    grid <- gds$findGridByShortName(shortName)
    if (is.null(grid)) {
        stop("Variable requested not found. Check variable nomenclature")
    }
    latLon <- getLatLonDomain(grid, lonLim, latLim)
    timePars <- getTimeDomain(grid, season, years, verifTime)
    levelPars <- getVerticalLevelPars(grid, level)
    mdArray <- makeSubset(grid, timePars$tRanges, levelPars$zRange, latLon)
    gds$close()
    if (!is.null(dictionary)) {
        ltb <- as.difftime(dic$lower_time_bound, format = "%H", units = "hours")
        utb <- as.difftime(dic$upper_time_bound, format = "%H", units = "hours")
        dateSliceStart <- as.POSIXlt(timePars$dateSlice - ltb)
        dateSliceEnd <- as.POSIXlt(timePars$dateSlice + utb)
        mdArray <- dictionaryTransform(dic, grid, timePars, mdArray) 
    } else {
        dateSliceStart <- timePars$dateSlice
        dateSliceEnd <- timePars$dateSlice + timePars$timeResInSeconds
    }
    message(paste("[",Sys.time(),"] - Done.", sep=""))
    return(list("VarName" = shortName, "isStandard" = isStandard, "Level" = levelPars$level, "Dates" = list("Start" = dateSliceStart,"End" = dateSliceEnd), "xyCoords" = latLon$xyCoords, "Data" = mdArray))
}
# End

 


