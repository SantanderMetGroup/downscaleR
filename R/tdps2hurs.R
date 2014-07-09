#' Conversion of dew point to relative humidity
#' 
#' @param tas vector of surface temperature (K)
#' @param tdps vector of surface dew point temperature (K)
#' @return vector of relative humidity, in %
#' @author J Bedia \email{joaquin.bedia@@gmail.com}, borrowing MatLab code from S. Herrera

tdps2hurs <- function(tas, tdps) {
      lv <- 2.5e+06
      Rv <- 461.5
      hurs <- 100 * exp((lv / Rv) * ((1 / tas) - (1 / tdps)))
      return(hurs)
}
# End 