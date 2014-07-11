


loadObservations.NetCDF <- function(dataset, var, lonLim, latLim) {
    ft <- J("ucar.nc2.constants.FeatureType")$valueOf("ANY_POINT") # Para sacar "ANY_POINT"
    nc <- J("ucar.nc2.dataset.NetcdfDataset")$openDataset(dataset)
    # Calculo de la longitud del eje tiempo 
    dims <- nc$getDimensions()
    dimNames <- unlist(strsplit(gsub("\\[|]|\\s|;", "", dims$toString()), split=","))
    grep("^time", dimNames)
    tsLen <- dims$get(as.integer(grep("^time", dimNames) - 1))$getLength()
    pds <- J("ucar.nc2.ft.FeatureDatasetFactoryManager")$wrap(ft, nc, NULL, .jnew("java.util.Formatter"))
    fc.list <- as.list(pds$getPointFeatureCollectionList())
    if (.jinstanceof(fc.list[[1]], "ucar.nc2.ft.StationTimeSeriesFeatureCollection")) {
	    stations.list <- as.list(fc.list[[1]]$getStations())
	    lons <- rep(NA, length(stations.list))
	    lats <- lons
	    alts <- lons
	    ids <- lons
	    for (i in 1:length(stations.list)) {
		    st <- fc.list[[1]]$getStationFeature(stations.list[[i]])
		    ids[i] <- st$getName() # Station code
		    lons[i] <- st$getLongitude() # Station lon
		    lats[i] <- st$getLatitude() # Station lat
		    alts[i] <- st$getAltitude() # Station alt
	    }
# 		st$getCalendarDateRange()$toString()
	########################
#	lonLim = c(-3,-2)
#	latLim = c(42,43)
	#######################
	if (is.null(lonLim) == FALSE) {
		latLon <- getLatLonDomainStations(lonLim, latLim, lons, lats)
		for (i in 1:length(latLon)) {
			assign(names(latLon[i]), latLon[[i]])
		}
	} else {
		stInd <- 1:length(stations.list)
	}
	######################
#	ids
#	stationID = c("000229")
	######################
	if (is.null(stationID) == FALSE) {
		stInd <- match(stationID, ids)
		if (any(is.na(stInd))) {
			stop("Non-existent stationID values found. Check data inventory")
		}
	}
	# Data retrieval 
	stDataList <- list()
	    for (i in 1:length(stInd)) {
#				tic(name = "total")
		    message(paste(Sys.time(), "Retrieving data from Station", ids[stInd[i]], "..."))
		    st <- fc.list[[1]]$getStationFeature(stations.list[[stInd[i]]])
		    tsDates <- rep(NA, tsLen)
		    tsValues <- tsDates
		    for (j in 1:tsLen) {
#					  tic(name = "tsLen")
			    st$hasNext()                 
			    pf <- J(st, "next")
#                 tsDates[j] <- pf$getObservationTimeAsCalendarDate()$toString()
			    tsValues[j] <- pf$getData()$getArray(var)$getStorage()
#						toc(name = "tsLen")
		    }      
		
#				  st$resetIteration()
#                   tsDates <- strptime(tsDates, format = "%Y-%m-%dT%H:%M:%SZ")
#                   stDataList[[length(stDataList) + 1]] <- list(tsDates, tsValues)
		    stDataList[[length(stDataList) + 1]] <- tsValues
#				  toc(name = "total")
	    }
	    str(stDataList)
	    stData <- do.call("cbind", stDataList)                  
    }
}


#             tsDates <- c()
#             tsValues <- c()
#             while(st$hasNext()) {
#                   pf <- J(st, "next")
#                   # Fecha
#                   tsDates <- pf$getObservationTimeAsCalendarDate()$toString())            
#                   # Valor
#                   tsValues <- c(tsValues, pf$getData()$getArray(var)$getStorage())
#             }   
#             st$resetIteration()
#             tsDates <- strptime(tsDates, format = "%Y-%m-%dT%H:%M:%SZ")
#             stDataList[[length(stDataList) + 1]] <- list(tsDates, tsValues)
#             }      
#             str(stDataList)

#### JAVA       
#nc = openDataset(location);
#PointDatasetImpl pointsDs = (PointDatasetImpl) FeatureDatasetFactoryManager.wrap(FeatureType.ANY_POINT, nc, null, new Formatter());
#for(FeatureCollection fc : pointsDs.getPointFeatureCollectionList()){
#	if(fc instanceof StationTimeSeriesFeatureCollection){
#		StationTimeSeriesFeatureCollection sfc = (StationTimeSeriesFeatureCollection) fc;
#		while(sfc.hasNext()){
#			// estacion
#			StationTimeSeriesFeature ts = sfc.next();
#			while(ts.hasNext()){
#				// elemento serie de la estacion
#				PointFeature feat = ts.next();
#				DateTime time = new DateTime(feat.getObservationTimeAsDate()).withZone(DateTimeZone.UTC);
#				if(time.equals(date)){
#					VariableDataBean vd = new VariableDataBean();
#					vd.setStation(ts);
#					vd.setValues(feat.getData().getJavaArrayFloat(layer.getId()));
#					vd.setDate(time);
#					data.put(ts, vd);
#					break;
#				}
#			}
#		}
#	}
#} 
#





