#' @title Grid containing NCEP reanalysis data of specific humidity at 850mb for the Iberian Peninsula.
#' @description The data are daily means, wintertime (DJF) period 1991-2010. 
#' @format A grid
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @name iberia_ncep_hus850
#' @examples
#' data(iberia_ncep_hus850)
#' plotMeanGrid(iberia_ncep_hus850)
NULL

#' @title Grid containing NCEP reanalysis data of sea-level pressure for the Iberian Peninsula.
#' @description The data correspond to the wintertime (DJF) period 1991-2010, and it consists of daily mean data, computed as the
#' mean of the four 6-hourly model outputs. 
#' @format A grid
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @name iberia_ncep_psl
#' @examples
#' data(iberia_ncep_psl)
#' plotMeanGrid(iberia_ncep_psl)
NULL

#' @title Grid containing NCEP reanalysis data of air temperature at 850mb for the Iberian Peninsula.
#' @description The data correspond to the wintertime (DJF) period 1991-2010, and it consists of daily mean data 
#' @format A grid
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @name iberia_ncep_ta850
#' @examples
#' data(iberia_ncep_ta850)
#' plotMeanGrid(iberia_ncep_ta850)
NULL

#' @title Multimember grid containing a seasonal forecast of maximum surface temperature for Europe
#' @description CFSv2 forecast of maximum daily temperature for July 2001 over Europe. Lead-month 2, first 9 members.
#' @format A multimember grid
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @references 
#' Saha, S. \emph{et al.}, 2013. The NCEP Climate Forecast System Version 2. J Clim 130925135638001. doi:10.1175/JCLI-D-12-00823.1
#' @name tasmax_forecast
#' @examples
#' data(tasmax_forecast)
#' plotMeanGrid(tasmax_forecast, multi.member = TRUE)
NULL

#' @title Multimember grid containing a seasonal forecast of minimum surface temperature for Europe
#' @description CFSv2 forecast of minimum daily temperature for July 2001 over Europe. Lead-month 2, first 9 members.
#' @format A multimember grid
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @references 
#' Saha, S. \emph{et al.}, 2013. The NCEP Climate Forecast System Version 2. J Clim 130925135638001. doi:10.1175/JCLI-D-12-00823.1
#' @name tasmin_forecast
#' @examples
#' data(tasmin_forecast)
#' plotMeanGrid(tasmin_forecast, multi.member = TRUE)
NULL

#' @title Multimember grid containing a seasonal forecast of precipitation for Europe
#' @description CFSv2 forecast of daily accumulated precipitation for July 2001 over Europe. Lead-month 2, first 9 members.
#' @format A multimember grid
#' @source \url{http://www.meteo.unican.es/ecoms-udg}
#' @references 
#' Saha, S. \emph{et al.}, 2013. The NCEP Climate Forecast System Version 2. J Clim 130925135638001. doi:10.1175/JCLI-D-12-00823.1
#' @name tp_forecast
#' @examples
#' data(tp_forecast)
#' plotMeanGrid(tp_forecast, multi.member = TRUE)
NULL

#' @title Grid containing NCEP reanalysis data of precipitation for the Iberian Peninsula.
#' @description NCEP_Iberia_tp is a grid object returned by loadECOMS from package loadeR.ECOMS
#' @format A grid
#' @name NCEP_Iberia_tp
#' @docType data
#' @keywords NCEP reanalysis
#' @source  subset of NCEP reanalysis data, which is accessible through the \strong{ECOMS User Data Gateway (ECOMS-UDG)} 
NULL

#' @title Grid containing NCEP reanalysis data of temperature for the Iberian Peninsula.
#' @description NCEP_Iberia_tas is a grid object returned by loadECOMS from package loadeR.ECOMS
#' @format A grid
#' @name NCEP_Iberia_tas
#' @docType data
#' @keywords NCEP reanalysis
#' @source  subset of NCEP reanalysis data, which is accessible through the \strong{ECOMS User Data Gateway (ECOMS-UDG)} 
NULL

#' @title Station daily precipitation dataset
#' @description Station data from the VALUE_ECA_86_v2 dataset containing daily precipitation for 82 stations in Europe.
#' @format Station data
#' @name VALUE_tp
#' @docType data
#' @source  Subset of VALUE station data. Full dataset is accessible 
#' for download in \strong{"http://meteo.unican.es/work/downscaler/data/VALUE_ECA_86_v2.tar.gz"}.
NULL

#' @title Station mean temperature dataset
#' @description Station data from the VALUE_ECA_86_v2 dataset containing daily mean temperature for 82 stations in Europe.
#' @format Station data
#' @name VALUE_tas
#' @docType data
#' @source  Subset of VALUE station data. Full dataset is accessible 
#' for download in \strong{"http://meteo.unican.es/work/downscaler/data/VALUE_ECA_86_v2.tar.gz"}.
NULL


#' @title Station daily precipitation data
#' @description Station data from the VALUE_ECA_86_v2 dataset containing daily precipitation for the Igueldo-SanSebastian station. 
#' @format Station data
#' @name VALUE_Igueldo_tp
#' @docType data
#' @source  Subset of VALUE station data. Full dataset is accessible 
#' for download in \strong{"http://meteo.unican.es/work/downscaler/data/VALUE_ECA_86_v2.tar.gz"}.
NULL

#' @title Station mean temperature data
#' @description Station data from the VALUE_ECA_86_v2 dataset containing daily mean temperature for the Igueldo-SanSebastian station. 
#' @format Station data
#' @name VALUE_Igueldo_tas
#' @docType data
#' @keywords VALUE station temperature Igueldo
#' @source  Subset of VALUE station data. Full dataset is accessible 
#' for download in \strong{"http://meteo.unican.es/work/downscaler/data/VALUE_ECA_86_v2.tar.gz"}.
NULL

#' @title Decada forecast example grid
#' @description Grid of decadal forecast data from the SPECS_GFDL_decadal dataset containing monthly temperature for the Iberian Peninsula. 
#' @format Grid data
#' @name tas_decadalForecast
#' @docType data
#' @keywords SPECS temperature decadal forecast
#' @source  Subset of the SPECS_GFDL_decadal dataset. This data is accessible 
#' for loading in the \strong{The User Data Gateway (UDG)}, which is the one stop shop for 
#' climate data access maintained by the Santander MetGroup. 
#' The UDG builds on the THREDDS Access Portal (UDG-TAP) which is the entry point for 
#' authentication and data access (more info in \link{https://meteo.unican.es/trac/wiki/udg}). 
NULL
