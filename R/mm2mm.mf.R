#' @title Constructor of a multimember multifield from various multimember fields
#' 
#' @description Constructs a multimember multifield from different multimember fields.
#'  Useful for the creation of multifields from the multimember fields returned by
#'  ecomsUDG.Raccess::loadECOMS in seasonal forecast downscaling applications.
#' 
#' @param mm.list Multimember list. A list containing the different input multimember fields.
#' 
#' @return A multimember multifield object
#' 
#' @details The function makes a number of checks in order to test the spatiotemporal compatibility of the input multi-member fields.
#'  Regarding the temporal concordance, it is implicitly assumeds that all temporal data from the different
#'  multimember fields correspond to the same time zone ("GMT"). The time zone itself is not important, as long as it is the
#'  same across datasets, because temporal consistency is checked on a daily basis, allowing the inclusion of predictors with
#'  different verification times and temporal aggregations. For instance, instantaneous geopotential at 12:00 is compatible with
#'  mean daily surface temperature, always that both variables correspond to the same days. Different time resolutions are not
#'  compatible and will return an error (for instance, 6-hourly data is incompatible with daily values, because their 
#'  respective time series have different lengths).
#' 
#' @export
#' 
#' @importFrom abind abind
#' 
#' @author J. bedia \email{joaquin.bedia@@gmail.com}
#' 
#' @seealso \code{\link[ecomsUDG.Raccess]{loadECOMS}}, \code{\link{loadMultifield}}
#' 


mm2mm.mf <- function(mm.list) {
      if (length(mm.list) < 2) {
            stop("The input must be a list of at least two multimember fields")
      }
      for (i in 2:length(mm.list)) {
            # Spatial test
            if (!identical(mm.list[[1]]$xyCoords, mm.list[[i]]$xyCoords)) {
                  stop("Input data is not spatially consistent")
            }
            # temporal test
            if (!identical(as.POSIXlt(mm.list[[1]]$Dates$start)$yday, as.POSIXlt(mm.list[[i]]$Dates$start)$yday) | !identical(as.POSIXlt(mm.list[[1]]$Dates$start)$year, as.POSIXlt(mm.list[[i]]$Dates$start)$year)) {
                  stop("Input data is not temporally consistent")
            }
            # data dimensionality
            if (!identical(dim(mm.list[[1]]$Data), dim(mm.list[[i]]$Data))) {
                  stop("Incompatible data array dimensions")
            }
            if (!identical(attr(mm.list[[1]]$Data, "dimensions"), attr(mm.list[[i]]$Data, "dimensions"))) {
                  stop("Inconsistent 'dimensions' attribute")
            }
      }
      aux.list <- lapply(1:length(mm.list), function(x) {mm.list[[x]]$Variable})
      varName <- unlist(lapply(1:length(aux.list), function(x) aux.list[[x]]$varName))      
      isStandard <- unlist(lapply(1:length(aux.list), function(x) aux.list[[x]]$isStandard))      
      level <- unlist(lapply(1:length(aux.list), function(x) ifelse(is.null(aux.list[[x]]$level), NA, aux.list[[x]]$level)))      
      aux.list <- NULL
      mm.list[[1]]$Variable <- list("varName" = varName, "isStandard" = isStandard, "level" = level)      
      mm.list[[1]]$Dates <- lapply(1:length(mm.list), function(x) {mm.list[[x]]$Dates})
      dimNames <- attr(mm.list[[1]]$Data, "dimensions")
      dimNames
      mm.list[[1]]$Data <- unname(do.call("abind", c(lapply(1:length(mm.list), function(x) {mm.list[[x]]$Data}), along = -1))) 
      attr(mm.list[[1]]$Data, "dimensions") <- c("var", dimNames)
      return(mm.list[[1]])
}
# End

# load("ignore/datasets/S4/wss12.Rdata")
# load("ignore/datasets/S4/hurs12.Rdata")
# load("ignore/datasets/S4/tas12.Rdata")
# load("ignore//datasets/S4//tpDA.Rdata", v = T)
# print(object.size(tas12.s4.15), units = "Mb")
# mmmf <- mm2mm.mf(mm.list = list(tas12.s4.15, hurs12.s4.15, tpDA.s4.15))
# 
# load("ignore//datasets/CFSv2//tasmax.Rdata", v = T)
# load("ignore//datasets/CFSv2//tasmin.Rdata", v = T)
# load("ignore//datasets/CFSv2//tp.Rdata", v = T)
# mmmf <- mm2mm.mf(list(tasmax.cfs.10, tasmin.cfs.9, tp.cfs.9))
# str(mmmf)

#¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬


