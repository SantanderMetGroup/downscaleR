# Description: Function to perform the spatial selection of stations based on lon-lat window definitions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getLatLonDomainStations <- function(lonLim, latLim, lons, lats) {
      lonLim <- lonLim
      latLim <- latLim
      lons <- lons
      lats <- lats
      if (length(lonLim) > 2 | length(latLim) > 2) {
            stop("Invalid definition of geographical position")
      }
      if (length(lonLim) != length(latLim)) {
            stop("Invalid definition of geographical position")
      }
      if (length(lonLim) == 2) {
            lonInd <- which(lons >= lonLim[1] & lons <= lonLim[2])      
            latInd <- which(lats >= latLim[1] & lats <= latLim[2])
            stInd <- intersect(lonInd, latInd)
      } else {
            stInd <- which.min(sqrt((lons-lonLim)^2 + (lats-latLim)^2))
      }
      return(list("stInd" = stInd, "stCoords" = cbind(lons, lats)[stInd, ]))
}
# End