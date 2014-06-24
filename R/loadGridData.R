loadGridData <- function(dataset, var, dictionary = TRUE, lonLim = NULL,
                         latLim = NULL, season = NULL, years = NULL, time = "none") {
    dataset <- dataset
    time <- match.arg(time, choices = c("none", "00", "06", "12", "18", "DD"))
    level <- findVerticalLevel(var)
    dic <- NULL
    if (isTRUE(dictionary)) {
        dicPath <- file.path(find.package("ecomsUDG.Raccess"), "dictionaries", paste(dataset,".dic", sep = ""))
        # devel
        # dicPath <- file.path("./inst/dictionaries", paste(dataset,".dic", sep = ""))
        dic <- dictionaryLookup(dicPath, var, time)
        shortName <- dic$short_name      
    } else {
        shortName <- var
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
    gds$close()
    message("[",Sys.time(),"]", " Done")
    return(out)
}      