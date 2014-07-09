#' Conversion of sea-level pressure to surface pressure
#' 
#' @param tas surface temperature
#' @param zs surface geopotential (m^2/s^2)
#' @param mslp sea-level pressure (Pa)
#' @return surface pressure (Pa)
#' @author Sixto Herrera
 
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