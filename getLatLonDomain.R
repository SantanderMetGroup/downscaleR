# ~~~~~~~~~~~~~~~ #
# getLatLonDomain
# ~~~~~~~~~~~~~~~ #
# Description
# ~~~~~~~~~~~~
# Low-level function used to retrieve the cartesian coordinates of a selected spatial slice -either single point locations or rectangular domains- from a grid coordinate system of the scientific data type "Grid" of netcdf-Java
# Arguments
# ~~~~~~~~~~~~
      # gridCoordinateSystem: coordinate system of the netcdf-Java class "Grid"
      # lonLim: argument passed by the loading function (see e.g. loadSeasonalForecast)
      # latLim: argument passed by the loading function (see e.g. loadSeasonalForecast)

# Notes
# ~~~~~~~~~~~~
# Implemented under netcdf-Java 4.3.18 / R 3.0.1/ rJava 0.9-4
# Juaco, 25 Sep 2013
# Modified 1 Oct 2013, dateline cross index introduced; netcdf-Java 4.3.18 / R 3.0.2/ rJava 0.9-4
# Modified 31 Oct 2013, account for bidimensional lonlat coordinate axes
# Modified 18 Jan 2014, bug fix, dateline crossing sub-routine relocated to be applied also in case lonLim is NULL
# Modified 23 Jan 2014, bug fix, selection fo the whole latitudinal domain for 2D lat definition, when latLim is NULL 
# ~~~~~~~~~~~~
getLatLonDomain <- function(gridCoordinateSystem, lonLim, latLim) {
      gcs <- gridCoordinateSystem
      if (length(lonLim) > 2 | length(latLim) > 2) {
    	      stop("Invalid coordinates. Check the spatial window definition")
      }
      lons <- gcs$getLonAxis()$getCoordValues()
      lons[which(lons > 180)] <- lons[which(lons > 180)] - 360
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	lonAxisShape <- gcs$getLonAxis()$getShape()
	if (length(lonAxisShape) == 2) {
		lons <- apply(t(matrix(lons, ncol = lonAxisShape[1])), 2, min)
	}
	#-------------------------------------------------------
	if (is.null(lonLim)) {
    	      lonSlice <- sort(lons)
    	      # lonInd <- (1:gcs$getLonAxis()$getShapeAll())[sort(lons, index.return = TRUE)$ix] - 1
	      lonInd <- (1:length(lons))[sort(lons, index.return = TRUE)$ix] - 1
    	      lonInd <- lonInd[order(lonSlice)]
    	      lonSlice <- sort(lonSlice)
      } else {
            if (length(lonLim) == 1) {
        	      lonInd <- which.min(abs(lons - lonLim)) - 1
        	      lonSlice <- lons[lonInd + 1]
        	      datelineCrossInd <- NULL
		}
            if (length(lonLim) == 2) {
			lonSlice <- subset(lons, lons >= lonLim[which.min(lonLim)] & lons <= lonLim[which.max(lonLim)])
                  lonInd <- which(lons %in% lonSlice) - 1
                  lonInd <- lonInd[order(lonSlice)]
                  lonSlice <- sort(lonSlice)
         	}
	}
      # dateline crossing 
      index <- c()
      for (i in 2:length(lonInd)) {
            index[i-1] <- lonInd[i] - lonInd[i-1]
      }
      datelineCrossInd <- which(index < 0)
      if (length(datelineCrossInd) == 0) {
            datelineCrossInd <- NULL
      }
      # Latitudes
      lats <- gcs$getLatAxis()$getCoordValues()
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	latAxisShape <- gcs$getLatAxis()$getShape()
	if (length(latAxisShape) == 2) {
	      lats <- apply(t(matrix(lats, ncol = lonAxisShape[1])), 1, min)
	}
	#-------------------------------------------------------
      if (is.null(latLim)) {
            latSlice <- lats
            # latInd <- (1:gcs$getLatAxis()$getShapeAll()) - 1
            latInd <- 1:length(lats) - 1
     	} else {
            if (length(latLim) == 1) {
                  latInd <- which.min(abs(lats - latLim)) - 1
                  latSlice <- lats[latInd + 1]
            }
            if (length(latLim) == 2) {
                  latSlice <- subset(lats, lats >= latLim[which.min(latLim)] & lats <= latLim[which.max(latLim)])
                  latInd <- which(lats %in% latSlice) - 1
            }
      }
	if (length(lonInd) < 1 | length(latInd) < 1) {
	      stop("Empty Spatial Domain.\nConsider expanding the boundaries or selecting a single point location")
	}
	gridPoints <- as.matrix(expand.grid("lon" = lonSlice, "lat" = latSlice))
	return(list("Grid" = gridPoints, "LatIndex" = latInd, "LonIndex" = lonInd, "DatelineCrossInd" = datelineCrossInd))
}
# End
