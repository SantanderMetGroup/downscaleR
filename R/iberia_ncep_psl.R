#' @title Field containing NCEP reanalysis data of sea-level pressure for the Iberian Peninsula.
#' 
#' @description The data correspond to the wintertime (DJF) period 1991-2010, and it consists of daily mean data, computed as the
#' mean of the four 6-hourly model outputs. 
#'
#' @format A field
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @name iberia_ncep_psl
#' @export
#' @seealso \code{\link{makeMultifield}}, \code{\link[ecomsUDG.Raccess]{loadECOMS}}
#' @examples
#' data(iberia_ncep_psl)
#' plotMeanField(iberia_ncep_psl)
NULL