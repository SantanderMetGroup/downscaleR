#' @title Load a field from a dataset
#' 
#' @description Load a user-defined spatio-temporal slice (a field) from a gridded dataset
#' 
#' @import rJava

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
#' \code{"00"}, \code{"06"}, \code{"12"} and \code{"18"}. If daily aggregated data are 
#' required use \code{"DD"}. If the requested variable is static (e.g. orography) it will be ignored. 
#' See details for time aggregation options
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
#' # Load air temperature at 850 millibar isobaric surface pressure level from the built-in NCEP dataset,
#' # for the Iberian Peninsula in summer (JJA):
#' ncep <- file.path(find.package("downscaleR"), "datasets/reanalysis/Iberia_NCEP/Iberia_NCEP.ncml")
#' field <- loadGridData(ncep, var = "ta@@85000", dictionary = TRUE, lonLim = c(-10,5),
#'    latLim = c(35.5, 44.5), season = 6:8, years = 1981:2010)
#' str(field)   
#' plotMeanField(field)
#' }
#' 

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
            stop("Variable '", shortName, "' not found\nCheck variable names using 'datasetInventory' and/or dictionary 'identifiers'")
      }
      latLon <- getLatLonDomain(grid, lonLim, latLim)
      out <- loadGridDataset(var, grid, dic, level, season, years, time, latLon)
      # Definition of projection
      proj <- grid$getCoordinateSystem()$getProjection()$toString()
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