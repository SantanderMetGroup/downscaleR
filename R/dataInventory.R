#' @name dataInventory
#' @description Function to provide a quick overview of a climate dataset
#'  (either stations or gridded data)
#' @param dataset A character string poiting to the target. Either a directory containing the data
#'  in the case of station data in standard ASCII format (see \code{\link{loadStationData}}) for details),
#'  or a target file (a NcML) in the case of other types of gridded data (reanalysis, gridded observations ...,
#'  see \code{\link{loadGridData}} for details).
#' @param return.stats Optional logical flag indicating if summary statistics of the dataset
#'  should be returned with the inventory. Only used for station data.
#'  
#' @return A list of components describing the variables and other characteristics of the target dataset.
#' @export
#' @author J Bedia \email{joaquin.bedia@@gmail.com}

dataInventory <- function(dataset, return.stats = FALSE) {
      rs <- return.stats
      message(paste("[", Sys.time(), "] Doing inventory ...", sep = ""))
      if (isTRUE(file.info(dataset)$isdir)) {
            out <- dataInventory.ASCII(dataset, rs)
      } else {
            out <- dataInventory.NetCDF(dataset)
      }
      message(paste("[", Sys.time(), "] Done.", sep = ""))
      return(out)
}
# End

