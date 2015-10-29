#' @title Load a field from a dataset
#' 
#' @description Load a user-defined spatio-temporal slice (a field) from a gridded dataset
#' 
#' @import rJava
#' 
#' @template templateParams
#' @param dictionary Default to TRUE, indicating that a dictionary is used and the .dic file is stored in the same path than the
#' dataset. If the .dic file is stored elsewhere, then the argument is the full path to the .dic file (including the extension,
#' e.g.: \code{"/path/to/the/dictionary_file.dic"}). This is the case for instance when the dataset is stored in a remote URL,
#' and we have a locally stored dictionary for that particular dataset. If FALSE no variable homogenization takes place,
#' and the raw variable, as originally stored in the dataset, will be returned. See details for dictionary specification.
#' @param time A character vector indicating the temporal filtering/aggregation 
#' of the output data. Default to \code{"none"}, which returns the original time 
#' series as stored in the dataset. For sub-daily variables, instantantaneous data at 
#' selected verification times can be filtered using one of the character strings 
#' \code{"00"}, \code{"03"}, \code{"06"}, \code{"09"}, \code{"12"}, \code{"15"},
#'  \code{"18"}, \code{"21"},and \code{"00"} when applicable. If daily aggregated data are 
#' required use \code{"DD"}. If the requested variable is static (e.g. orography) it will be ignored. 
#' See the next arguments for time aggregation options.
#' @param aggr.d Character string. Function of aggregation of sub-daily data for daily data calculation. 
#' Currently accepted values are \code{"none"}, \code{"mean"}, \code{"min"}, \code{"max"} and \code{"sum"}.
#' @param aggr.m Same as \code{aggr.d}, bun indicating the aggregation function to compute monthly from daily data.
#' If \code{aggr.m = "none"} (the default), no monthly aggregation is undertaken.
#' 
#' @template templateReturnGridData
#' @template templateDicDetails  
#' @template templateGeolocation
#' @export
#' 
#' @note The possible values of argument \code{dictionary} are slightly different from those from the
#'  \code{\link[ecomsUDG.Raccess]{loadECOMS}} function in package \pkg{ecomsUDG.Raccess}, basically because
#'  the latter does not accept user-defined dictionaries at arbitrary file paths.
#'    
#' @author J. Bedia \email{joaquin.bedia@@gmail.com}
#' @family loading
#' @family loading.grid
#' @family homogenization
#' 
#' @examples \dontrun{
#' #Download dataset
#' dir.create("mydirectory")
#' download.file("http://meteo.unican.es/work/downscaler/data/Iberia_NCEP.tar.gz", 
#' destfile = "mydirectory/Iberia_NCEP.tar.gz")
#' # Extract files from the tar.gz file
#' untar("mydirectory/NCEP_Iberia.tar.gz", exdir = "mydirectory")
#' # First, the path to the ncml file is defined:
#' ncep <- "mydirectory/Iberia_NCEP/Iberia_NCEP.ncml"
#' # Load air temperature at 850 millibar isobaric surface pressure level from the built-in
#' # NCEP dataset, for the Iberian Peninsula in summer (JJA):
#' field <- loadGridData(ncep, var = "ta@@850", dictionary = TRUE, lonLim = c(-10,5),
#'    latLim = c(35.5, 44.5), season = 6:8, years = 1981:2010)
#' str(field)   
#' plotMeanField(field)
#' # Calculation of monthly mean temperature:
#' field.mm <- loadGridData(ncep, var = "ta@@850", dictionary = TRUE, lonLim = c(-10,5),
#'                          latLim = c(35.5, 44.5), season = 6:8,
#'                          years = 1981:2010, aggr.m = "mean")
#' str(field.mm)
#' 
#' # Same but using the original variable (not homogenized via dictionary):
#' di <- dataInventory(ncep)
#' names(di)
#' data(vocabulary)
#' vocabulary
#' # Variable is named 'T', instead of the standard name 'ta' in the vocabulary
#' # Vertical level is indicated using the '@@' symbol:
#' non.standard.field <- loadGridData(ncep, var = "T@@850", dictionary = FALSE, lonLim = c(-10,5),
#'                                   latLim = c(35.5, 44.5), season = 6:8, 
#'                                   years = 1981:2010, aggr.m = "mean")
#' str(non.standard.field$Variable)
#' # Note the units are now in Kelvin, as originally stored
#' plotMeanField(non.standard.field)
#' }
#' 

loadGridData <- function(dataset, var, dictionary = TRUE, lonLim = NULL,
                         latLim = NULL, season = NULL, years = NULL, time = "none",
                         aggr.d = "none", aggr.m = "none") {
      time <- match.arg(time, choices = c("none","00","03","06","09","12","15","18","21","DD"))
      aggr.d <- match.arg(aggr.d, choices = c("none", "mean", "min", "max", "sum"))
      if (time != "DD" & aggr.d != "none") {
            aggr.d <- "none"
            message("NOTE: Argument 'aggr.d' ignored as 'time' was set to ", time)
      }
      aggr.m <- match.arg(aggr.m, choices = c("none", "mean", "min", "max", "sum"))
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
      if (!is.null(season) & (min(season) < 1 | max(season) > 12)) {
            stop("Invalid season definition")
      }
      gds <- J("ucar.nc2.dt.grid.GridDataset")$open(dataset)
      grid <- gds$findGridByShortName(shortName)
      if (is.null(grid)) {
            stop("Variable '", shortName, "' not found\nCheck variable names using 'datasetInventory' and/or dictionary 'identifiers'")
      }
      latLon <- getLatLonDomain(grid, lonLim, latLim)
      proj <- grid$getCoordinateSystem()$getProjection()
      if (!proj$isLatLon()){
        nc <- gds$getNetcdfDataset()
        lonAxis <- nc$findVariable('lon')
        auxLon <- t(matrix(data = lonAxis$getCoordValues(), nrow = lonAxis$getShape()[2], ncol = lonAxis$getShape()[1]))
        latAxis <- nc$findVariable('lat')
        auxLat <- t(matrix(data = latAxis$getCoordValues(), nrow = latAxis$getShape()[2], ncol = latAxis$getShape()[1]))
        if (is.null(lonLim)){
          lonLim <- c(min(auxLon),max(auxLon))
        }
        if (is.null(latLim)){
          latLim <- c(min(auxLat),max(auxLat))
        }
        if (length(lonLim) == 1 | length(latLim) == 1) {
          ind.x <- which.min(abs(auxLon - lonLim))
          ind.y <- which.min(abs(auxLat - latLim))
          pointXYindex <- c(ind.y,ind.x)
          latLon$xyCoords$x <- grid$getCoordinateSystem()$getXHorizAxis()$getCoordValues()[ind.x]
          latLon$xyCoords$y <- grid$getCoordinateSystem()$getYHorizAxis()$getCoordValues()[ind.y]
          latLon$xyCoords$lon <- auxLon[ind.y,ind.x]
          latLon$xyCoords$lat <- auxLat[ind.y,ind.x]
        }else{
          auxDis <- sqrt((auxLon - lonLim[1])^2+(auxLat - latLim[1])^2)
          llrowCol <- arrayInd(which.min(auxDis), dim(auxDis))
          auxDis <- sqrt((auxLon - lonLim[2])^2+(auxLat - latLim[2])^2)
          urrowCol <- arrayInd(which.min(auxDis), dim(auxDis))
          auxDis <- sqrt((auxLon - lonLim[1])^2+(auxLat - latLim[2])^2)
          ulrowCol <- arrayInd(which.min(auxDis), dim(auxDis))
          auxDis <- sqrt((auxLon - lonLim[2])^2+(auxLat - latLim[1])^2)
          lrrowCol <- arrayInd(which.min(auxDis), dim(auxDis))
          llrowCol <- c(min(c(llrowCol[1],lrrowCol[1])), min(c(llrowCol[2],ulrowCol[2])))
          urrowCol <- c(max(c(ulrowCol[1],urrowCol[1])),max(c(lrrowCol[2],ulrowCol[2])))
          latLon$xyCoords$x <- grid$getCoordinateSystem()$getXHorizAxis()$getCoordValues()[llrowCol[2]:urrowCol[2]]
          latLon$xyCoords$y <- grid$getCoordinateSystem()$getYHorizAxis()$getCoordValues()[llrowCol[1]:urrowCol[1]]
          latLon$xyCoords$lon <- auxLon[llrowCol[1]:urrowCol[1],llrowCol[2]:urrowCol[2]]
          latLon$xyCoords$lat <- auxLat[llrowCol[1]:urrowCol[1],llrowCol[2]:urrowCol[2]]
        }
        latLon$lonRanges <- .jnew("ucar/ma2/Range", as.integer(llrowCol[2]-1), as.integer(urrowCol[2]-1))
        latLon$latRanges <- .jnew("ucar/ma2/Range", as.integer(llrowCol[1]-1), as.integer(urrowCol[1]-1))
      }
      out <- loadGridDataset(var, grid, dic, level, season, years, time, latLon, aggr.d, aggr.m)
      # Definition of projection
      proj <- proj$toString()
      attr(out$xyCoords, which = "projection") <- proj
      # Dimension ordering
      tab <- c("time", "level", "lat", "lon")
      x <- attr(out$Data, "dimensions")
      if (length(x) > 1) {
        b <- na.exclude(match(tab, x))
        dimNames <- attr(out$Data, "dimensions")[b]
        out$Data <- aperm(out$Data, perm = b)    
        attr(out$Data, "dimensions")  <- dimNames
      }
      # Source Dataset and other metadata 
      attr(out, "dataset") <- dataset
      gds$close()
      message("[",Sys.time(),"]", " Done")
      return(out)
}     
# End