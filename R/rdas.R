#' @title Field containing NCEP reanalysis data of specific humidity at 850mb for the Iberian Peninsula.
#' 
#' @description The data correspond to the wintertime (DJF) period 1991-2010, and it consists of instantaneous data at 12:00 UTC 
#'
#' @format A field
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @name iberia_ncep_hus850
#' @seealso \code{\link{makeMultiField}}, \code{\link[ecomsUDG.Raccess]{loadECOMS}}
#' @examples
#' data(iberia_ncep_hus850)
#' plotMeanField(iberia_ncep_hus850)
NULL

#' @title Field containing NCEP reanalysis data of sea-level pressure for the Iberian Peninsula.
#' 
#' @description The data correspond to the wintertime (DJF) period 1991-2010, and it consists of daily mean data, computed as the
#' mean of the four 6-hourly model outputs. 
#'
#' @format A field
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @name iberia_ncep_psl
#' @seealso \code{\link{makeMultiField}}, \code{\link[ecomsUDG.Raccess]{loadECOMS}}
#' @examples
#' data(iberia_ncep_psl)
#' plotMeanField(iberia_ncep_psl)
NULL

#' @title Field containing NCEP reanalysis data of air temperature at 850mb for the Iberian Peninsula.
#' 
#' @description The data correspond to the wintertime (DJF) period 1991-2010, and it consists of instantaneous data at 12:00 UTC 
#'
#' @format A field
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @name iberia_ncep_ta850
#' @seealso \code{\link{makeMultiField}}, \code{\link[ecomsUDG.Raccess]{loadECOMS}}
#' @examples
#' data(iberia_ncep_ta850)
#' plotMeanField(iberia_ncep_ta850)
NULL

#' @title Multimember field containing a seasonal forecast of maximum surface temperature for Europe
#' 
#' @description CFSv2 forecast of maximum daily temperature for July 2001 over Europe. Lead-month 2, first 9 members.
#'
#' @format A multimember field
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @references 
#' Saha, S. \emph{et al.}, 2013. The NCEP Climate Forecast System Version 2. J Clim 130925135638001. doi:10.1175/JCLI-D-12-00823.1
#' @name tasmax_forecast
#' @seealso \code{\link{makeMultiField}}, \code{\link[ecomsUDG.Raccess]{loadECOMS}}
#' @examples
#' data(tasmax_forecast)
#' plotMeanField(tasmax_forecast, multi.member = TRUE)
NULL

#' @title Multimember field containing a seasonal forecast of minimum surface temperature for Europe
#' 
#' @description CFSv2 forecast of minimum daily temperature for July 2001 over Europe. Lead-month 2, first 9 members.
#'
#' @format A multimember field
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @references 
#' Saha, S. \emph{et al.}, 2013. The NCEP Climate Forecast System Version 2. J Clim 130925135638001. doi:10.1175/JCLI-D-12-00823.1

#' @name tasmin_forecast
#' @seealso \code{\link{makeMultiField}}, \code{\link[ecomsUDG.Raccess]{loadECOMS}}
#' @examples
#' data(tasmin_forecast)
#' plotMeanField(tasmin_forecast, multi.member = TRUE)
NULL

#' @title Multimember field containing a seasonal forecast of precipitation for Europe
#' 
#' @description CFSv2 forecast of daily accumulated precipitation for July 2001 over Europe. Lead-month 2, first 9 members.
#'
#' @format A multimember field
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @references 
#' Saha, S. \emph{et al.}, 2013. The NCEP Climate Forecast System Version 2. J Clim 130925135638001. doi:10.1175/JCLI-D-12-00823.1

#' @name tp_forecast
#' @seealso \code{\link{makeMultiField}}, \code{\link[ecomsUDG.Raccess]{loadECOMS}}
#' @examples
#' data(tp_forecast)
#' plotMeanField(tp_forecast, multi.member = TRUE)
NULL
