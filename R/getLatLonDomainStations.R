getLatLonDomainStations <- function(lonLim, latLim, lons, lats) {
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
            message("[", Sys.time(),"] Closest station located at ", round(min(sqrt((lons-lonLim)^2 + (lats-latLim)^2)), digits=4), " spatial units from the specified [lonLim,latLim] coordinate") 
      
      }
      return(list("stInd" = stInd, "stCoords" = as.matrix(cbind(lons[stInd], lats[stInd]))))
}
# End