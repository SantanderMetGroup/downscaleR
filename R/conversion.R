#' @title Specific humidity from relative humidity
#' @description Derive specific humidity from relative humidity
#' 
#' @param tas surface air temperature (K)
#' @param ps surface pressure (Pa)
#' @param hurs surface relative humidity (\%)
#' @return Specific humidity (kg.kg-1)
#' @references Bohren & Albrecht (2000) Atmospheric thermodynamics. Oxford University Press. 402 pp
#' @author S. Herrera \email{sixto@@predictia.es}
#' @keywords internal
#' @export
#' @family conversion

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
#################################################################################

#' @title Sea-level pressure to surface pressure
#' @description Conversion of sea-level pressure to surface pressure
#' 
#' @param tas surface temperature
#' @param zs surface geopotential (m^2/s^2)
#' @param mslp sea-level pressure (Pa)
#' @return surface pressure (Pa)
#' @author S. Herrera \email{sixto@@predictia.es}
#' @keywords internal
#' @family conversion
#' @export

mslp2ps <- function(tas, zs, mslp) {
      Rd <- 287.058 # dry air constant J/(K kg)
      GammaST <- 0.0065 #(dT/dz)^st standard atmosphere vertical gradient of the temperature in the troposphere (0.0065 (K/m^-1))
      g <- 9.80665 # gravity (m/s^2)
      To <- tas + GammaST * zs / g
      ind <- which(abs(zs) >= 0.001)
      ind1 <- intersect(intersect(which(To > 290.5), which(tas <= 290.5)), ind)
      auxGamma <- mslp
      ps <- mslp
      auxGamma[ind1] <- g * (290.5 - tas[ind1]) / zs[ind1]
      ind <- setdiff(ind, ind1)
      ind1 <- intersect(intersect(which(To > 290.5), which(tas > 290.5)), ind)
      auxGamma[ind1] <- 0
      tas[ind1] <- 0.5 * (255 + tas[ind1])
      ind <- setdiff(ind, ind1)
      ind1 <- intersect(which(tas < 255), ind)
      auxGamma[ind1] <- GammaST
      tas[ind1] <- 0.5 * (255 + tas[ind1])
      ind <- setdiff(ind, ind1)
      auxGamma[ind] <- GammaST
      ind <- which(abs(zs) >= 0.001)
      ps[ind] <- mslp[ind] * exp((-zs[ind] / (Rd * tas[ind])) * (1 - 0.5 * (auxGamma[ind] * zs[ind]) / (g * tas[ind]) + (1 / 3) * ((auxGamma[ind] * zs[ind]) / (g * tas[ind])) ^ 2))
      tas <- NULL
      zs <- NULL
      mslp <- NULL
      auxGamma <- NULL
      return(ps)
}
# End
#################################################################################

#' @title Dew point to relative humidity
#' @description Conversion of dew point to relative humidity
#' 
#' @param tas vector of surface temperature (K)
#' @param tdps vector of surface dew point temperature (K)
#' @return vector of relative humidity (\%)
#' @author J Bedia \email{joaquin.bedia@@gmail.com}, borrowing MatLab code from S. Herrera
#' @keywords internal
#' @family conversion
#' @export

tdps2hurs <- function(tas, tdps) {
      lv <- 2.5e+06
      Rv <- 461.5
      hurs <- 100 * exp((lv / Rv) * ((1 / tas) - (1 / tdps)))
      return(hurs)
}
# End 