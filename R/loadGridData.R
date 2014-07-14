#' @title Load a gridded dataset
#' 
#' @description Load a user-defined spatio-temporal slice from a gridded dataset
#' 
#' @import rJava

#' @param dataset
#' @param var
#' @param dictionary
#' @param lonLim
#' @param latLim
#' @param season
#' @param years
#' @param time
#' 
#' @return a list of elements
#' 
#' @export
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @aliases loading

loadGridData <- function(dataset, var, dictionary = TRUE, lonLim = NULL,
                         latLim = NULL, season = NULL, years = NULL, time = "none") {
    time <- match.arg(time, choices = c("none", "00", "06", "12", "18", "DD"))
    aux.level <- findVerticalLevel(var)
    var <- aux.level$var
    level <- aux.level$level
    # Dictionary lookup
    if (dictionary == FALSE) {
        dic <- NULL
        shortName <- var
    } else {
        if (isTRUE(dictionary)) {
            dicPath <- gsub("ncml$", "dic", dataset)
        }
        if (is.character(dictionary)) {
            dicPath <- dictionary
        }
        dic <- dictionaryLookup(dicPath, var, time)
        shortName <- dic$short_name          
    }
    if (is.null(season)) {
        season <- 1:12
    }
    if (min(season) < 1 | max(season) > 12) {
        stop("Invalid season definition")
    }
    gds <- J("ucar.nc2.dt.grid.GridDataset")$open(dataset)
    grid <- gds$findGridByShortName(shortName)
    if (is.null(grid)) {
        stop("Variable requested not found")#.\nCheck variables using 'datasetInventory'")
    }
    latLon <- getLatLonDomain(grid, lonLim, latLim)
    out <- loadGridDataset(var, grid, dic, level, season, years, time, latLon)
    # Definition of projection
    proj <- grid$getCoordinateSystem()$getProjection()$toString()
    attr(out$xyCoords, which = "projection") <- proj
    gds$close()
    message("[",Sys.time(),"]", " Done")
    return(out)
}      