#' @title Constructor of a multifield from various fields, supporting multimember fields.
#' 
#' @description Constructs a (possibly multimember) multifield from different (multimember) fields.
#' 
#' @param ... Input fields to form the multifield. These must be compatible in time and space (see details).
#' @param spatial.tolerance numeric. Coordinate differences smaller than \code{spatial.tolerance} will be considered equal 
#' coordinates. Default to 0.001 --assuming that degrees are being used it seems a reasonable rounding error after interpolation--.
#' This value is passed to the \code{\link{identical}} function to check for spatial consistency of the input fields.
#' 
#' @return A (multimember) multifield object encompassing the different input (multimember) fields
#' 
#' @details The function makes a number of checks in order to test the spatiotemporal compatibility of the input multi-member fields.
#'  Regarding the temporal concordance, it is implicitly assumed that all temporal data from the different
#'  multimember fields correspond to the same time zone ("GMT"). The time zone itself is not important, as long as it is the
#'  same across datasets, because temporal consistency is checked on a daily basis (not hourly), allowing the inclusion of 
#'  predictors with different verification times and temporal aggregations. For instance, instantaneous geopotential at 12:00 
#'  is compatible with mean daily surface temperature, always that both variables correspond to the same days. Different time 
#'  resolutions are not compatible and will return an error (for instance, 6-hourly data is incompatible with daily values, 
#'  because their respective time series for a given season have different lengths).
#'  
#'  The spatial consistency of the input fields is also checked. In order to avoid possible errors from the user, the spatial
#'   consistency (i.e., equal XY coordinates) of the input fields must be ensured before attempting the creation of the multifield,
#'   otherwise giving an error. This can be achieved either through the specification of the same 'lonLim' and 'latLim' argument
#'   values when loading the fields, or using the \code{\link{interpGridData}} interpolator in conjuntion with the \code{\link{getGrid}}
#'   method.
#'  
#'  
#' @note Multifield can not be passed to the interpolator \code{\link{interpGridData}} directly. Instead, the 
#' multimember fields should be interpolated individually prior to multifield construction. 
#' 
#' @export
#' 
#' @importFrom abind abind
#' 
#' @author J. bedia \email{joaquin.bedia@@gmail.com}
#' 
#' @seealso \code{\link{loadGridData}} and \code{\link[ecomsUDG.Raccess]{loadECOMS}} for loading fields (the latter also for loading
#' multimember fields), \code{\link{loadMultiField}}, which directly loads a multifield. \code{\link{interpGridData}} and \code{\link{getGrid}}
#' for spatial consistency of input fields.
#' 
#' @examples
#' 
#' # Creation of a multifield from three different fields:
#' data(iberia_ncep_ta850)
#' data(iberia_ncep_hus850)
#' data(iberia_ncep_psl)
#' # An example of different temporal aggregations, temporally compatible: sea-level pressure is a daily mean,
#' # while specific humidity and air temperature (850 mb surface isobaric pressure level) are instantaneous data verifying at
#' # 12:00 UTC:
#' # air temperature
#' range(iberia_ncep_ta850$Dates$start)
#' range(iberia_ncep_ta850$Dates$end) # start and end are identical (instantaneous)
#' # sea-level pressure
#' range(iberia_ncep_psl$Dates$start)
#' range(iberia_ncep_psl$Dates$end) # start and end differ in 24 h (daily mean)
#' mf <- makeMultiField(iberia_ncep_hus850, iberia_ncep_psl, iberia_ncep_ta850)
#' # The new object inherits the global attributes from the first field, as it is assumed
#' # that all input fields come from the same data source:
#' attributes(mf)
#' # The data structure has now one additional dimension ("var"), along which the data arrays have been binded:
#' str(mf$Data)
#' plotMeanField(mf)
#' 
#' # Example of multimember multifield creation from several multimember fields:
#' # Load three different multimember fields with the same spatiotemporal ranges:
#' data(tasmax_forecast)
#' data(tasmin_forecast)
#' data(tp_forecast)
#' mm.mf <- makeMultiField(tasmax_forecast, tasmin_forecast, tp_forecast)
#' # 'plotMeanField' can just handle the multi-member mean for each variable in this case:
#' plotMeanField(mm.mf)
#' 

makeMultiField <- function(..., spatial.tolerance = 1e-3) {
      field.list <- list(...)
      if (length(field.list) < 2) {
            stop("The input must be a list of at least two multimember fields")
      }
      tol <- spatial.tolerance
      for (i in 2:length(field.list)) {
            # Spatial test
            if (!isTRUE(all.equal(field.list[[1]]$xyCoords, field.list[[i]]$xyCoords, check.attributes = FALSE, tolerance = tol))) {
                  stop("Input data is not spatially consistent")
            }
            # temporal test
            if (!identical(as.POSIXlt(field.list[[1]]$Dates$start)$yday, as.POSIXlt(field.list[[i]]$Dates$start)$yday) | !identical(as.POSIXlt(field.list[[1]]$Dates$start)$year, as.POSIXlt(field.list[[i]]$Dates$start)$year)) {
                  stop("Input data is not temporally consistent")
            }
            # data dimensionality
            if (!identical(dim(field.list[[1]]$Data), dim(field.list[[i]]$Data))) {
                  stop("Incompatible data array dimensions")
            }
            if (!identical(attr(field.list[[1]]$Data, "dimensions"), attr(field.list[[i]]$Data, "dimensions"))) {
                  stop("Inconsistent 'dimensions' attribute")
            }
      }
      aux.list <- lapply(1:length(field.list), function(x) {field.list[[x]]$Variable})
      varName <- unlist(lapply(1:length(aux.list), function(x) aux.list[[x]]$varName))      
      isStandard <- unlist(lapply(1:length(aux.list), function(x) aux.list[[x]]$isStandard))      
      level <- unlist(lapply(1:length(aux.list), function(x) ifelse(is.null(aux.list[[x]]$level), NA, aux.list[[x]]$level)))      
      aux.list <- NULL
      field.list[[1]]$Variable <- list("varName" = varName, "isStandard" = isStandard, "level" = level)      
      field.list[[1]]$Dates <- lapply(1:length(field.list), function(x) {field.list[[x]]$Dates})
      dimNames <- attr(field.list[[1]]$Data, "dimensions")
      field.list[[1]]$Data <- unname(do.call("abind", c(lapply(1:length(field.list), function(x) {field.list[[x]]$Data}), along = -1))) 
      attr(field.list[[1]]$Data, "dimensions") <- c("var", dimNames)
      return(field.list[[1]])
}
# End



