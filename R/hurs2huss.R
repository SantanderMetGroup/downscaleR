#' Derive specific humidity from relative humidity
#' 
#' @param tas surface air temperature (K)
#' @param ps surface pressure (Pa)
#' @param hurs surface relative humidity (%)
#' @return Specific humidity (kg.kg-1)
#' @references Bohren & Albrecht (2000) Atmospheric thermodynamics. Oxford University Press. 402 pp
#' @author Sixto Herrera

hurs2huss <- function(tas, ps, hurs) {
      Rd <- 287.058 # dry air constant J/(K kg)
      Rv <- 461.5 # J/(K kg)
      T0 <- 273.15 # (K)
      es0 <- 611
      es <- tas
      # Saturation pressure (Pa) over ice and water respectively (Bohren & Albrecht 2000, pp 197-200)
      iceMask <- which(tas < T0)
      es[iceMask] <- es0 * exp((6293 / T0) - (6293 / tas[iceMask]) - 0.555 * log(abs(tas[iceMask] / T0)))
      waterMask <- which(tas >= T0)
      es[waterMask] <- es0 * exp((6808 / T0) - (6808 / tas[waterMask]) - 5.09 * log(abs(tas[waterMask] / T0)))
      tas <- NULL
      ws <- (Rd / Rv) * (es / (ps - es))
      ps <- NULL
      es <- NULL
      w <- ws * hurs * 0.01
      ws <- NULL
      hurs <- NULL
      huss <- w / (1 + w)
      w <- NULL
      return(huss)
}
# End