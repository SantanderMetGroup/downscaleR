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

#' @title Field containing NCEP reanalysis data of precipitation for the Iberian Peninsula.
#' @description NCEP_Iberia_tp is a field object returned by loadECOMS 
#' function (package ecomsUDG.Raccess):
#' season = c(12,1,2),  
#' years = 1991-2000, 
#' var = "tp", 
#' lonLim = c(-10,5),  
#' latLim= c(34,44), 
#' time = "DD", 
#' aggr.d = "sum".
#' @format A field
#' @name NCEP_Iberia_tp
#' @docType data
#' @keywords NCEP reanalysis
#' @source  subset of NCEP reanalysis data, which is accessible 
#' through the \strong{ECOMS User Data Gateway (ECOMS-UDG)} 
#' @seealso \code{\link{makeMultiField}}, \code{\link[ecomsUDG.Raccess]{loadECOMS}}


NULL

#' @title Field containing NCEP reanalysis data of temperature for the Iberian Peninsula.
#' @description NCEP_Iberia_tas is a field object returned by loadECOMS 
#' function (package ecomsUDG.Raccess):
#' season = c(12,1,2),  
#' years = 1991-2000, 
#' var="tas", 
#' lonLim = c(-10,5),  
#' latLim= c(34,44), 
#' time = "DD", 
#' aggr.d = "min".
#' @format A field
#' @name NCEP_Iberia_tas
#' @docType data
#' @keywords NCEP reanalysis
#' @source  subset of NCEP reanalysis data, which is accessible 
#' through the \strong{ECOMS User Data Gateway (ECOMS-UDG)} 
#' @seealso \code{\link{makeMultiField}}, \code{\link[ecomsUDG.Raccess]{loadECOMS}}


NULL


#' @title Station data from the VALUE_ECA_86_v2 dataset containing daily precipitation for 82 stations in Europe.
#' @description 
#' season = c(12,1,2),  
#' years = 1991-2000, 
#' var= "precip".
#' @format Station data
#' @name VALUE_tp
#' @docType data
#' @keywords VALUE station precipitation 
#' @source  Subset of VALUE station data. Full dataset is accessible 
#' for download in \strong{"http://meteo.unican.es/work/downscaler/data/VALUE_ECA_86_v2.tar.gz"}.
#' @seealso \code{\link{NCEP_Iberia_tp}}, \code{\link{VALUE_Igueldo_tp}}


NULL


#' @title Station data from the VALUE_ECA_86_v2 dataset containing daily mean temperature for 82 stations in Europe.
#' @description 
#' season = c(12,1,2),  
#' years = 1991-2000, 
#' var= "tmean".
#' @format Station data
#' @name VALUE_tas
#' @docType data
#' @keywords VALUE station temperature
#' @source  Subset of VALUE station data. Full dataset is accessible 
#' for download in \strong{"http://meteo.unican.es/work/downscaler/data/VALUE_ECA_86_v2.tar.gz"}.
#' @seealso \code{\link{NCEP_Iberia_tas}}, \code{\link{VALUE_Igueldo_tas}}


NULL


#' @title Station data from the VALUE_ECA_86_v2 dataset containing daily precipitation for the Igueldo-SanSebastian station.
#' @description 
#' season = c(12,1,2),  
#' years = 1991-2000, 
#' stationID = "000234"
#' var= "precip".
#' @format Station data
#' @name VALUE_Igueldo_tp
#' @docType data
#' @keywords VALUE station precipitation Igueldo
#' @source  Subset of VALUE station data. Full dataset is accessible 
#' for download in \strong{"http://meteo.unican.es/work/downscaler/data/VALUE_ECA_86_v2.tar.gz"}.
#' @seealso \code{\link{NCEP_Iberia_tp}}, \code{\link{VALUE_tp}}


NULL

#' @title Station data from the VALUE_ECA_86_v2 dataset containing daily mean temperature for the Igueldo-SanSebastian station.
#' @description 
#' season = c(12,1,2),  
#' years = 1991-2000, 
#' stationID = "000234"
#' var= "tmean".
#' @format Station data
#' @name VALUE_Igueldo_tas
#' @docType data
#' @keywords VALUE station temperature Igueldo
#' @source  Subset of VALUE station data. Full dataset is accessible 
#' for download in \strong{"http://meteo.unican.es/work/downscaler/data/VALUE_ECA_86_v2.tar.gz"}.
#' @seealso \code{\link{NCEP_Iberia_tas}}, \code{\link{VALUE_tas}}


NULL
