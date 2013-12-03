dataInventory.NetCDF <- function(dataset) {
      dataset <- dataset
      # jString <- .jnew("java/lang/String", dataset)
      gds <- J("ucar.nc2.dt.grid.GridDataset")$open(dataset)
      varNames <- unlist(strsplit(gsub("\\[|]|\\s", "", gds$getGrids()$toString()), ","))
      if (length(varNames) == 0) {
            # datasetType = "PointDataset"
            ncds = gds$getNetcdfDataset()
            varList = .jevalArray(ncds$getVariables()$toArray())
            varNames = rep(NA, length(varList))
            for (k in 1:length(varList)) {
                  varNames[k] = varList[[k]]$getShortName()
            }
            ## Station IDs -----------------------------------------------------------------------
            stationIds = unlist(strsplit(varList[[grep("station", varNames)]]$read()$toString(),","))
            ## Coordinates -----------------------------------------------------------------------
            coordSysList = .jevalArray(ncds$getCoordinateAxes()$toArray())
            coordNames = rep(NA, length(coordSysList))
            for (k in 1:length(coordSysList)) {
                  coordNames[k] = coordSysList[[k]]$getShortName()
            }
            ## Longitude and Latitude ------------------------------------------------------------
            lonUnits = coordSysList[[grep("longitude", coordNames)]]$getUnitsString()
            lon = coordSysList[[grep("longitude", coordNames)]]$getCoordValues()
            latUnits = coordSysList[[grep("latitude", coordNames)]]$getUnitsString()
            lat = coordSysList[[grep("latitude", coordNames)]]$getCoordValues()
            lonLatList = list("Units" = paste(c(latUnits, lonUnits), collapse = ","), "Values" = cbind(lat,lon))
            ## Altitude --------------------------------------------------------------------------
            altUnits = coordSysList[[grep("altitude", coordNames)]]$getUnitsString()
            altValues = coordSysList[[grep("altitude", coordNames)]]$getCoordValues()
            altList = list("Units" = altUnits, "Values" = altValues)
            ## Times -----------------------------------------------------------------------------
            timeUnits = coordSysList[[grep("time", coordNames)]]$getUnitsString()
            timeStep = coordSysList[[grep("time", coordNames)]]$getIncrement()
            timeRange = c(coordSysList[[grep("time", coordNames)]]$getMinValue(), coordSysList[[grep("time", coordNames)]]$getMaxValue())
            timeList = list("Units" = timeUnits, "TimeStep" = timeStep, "TimeRange" = timeRange)
            ## Rest of variables -----------------------------------------------------------------
            climVars = varNames[which(varNames %in% c(coordNames, grep("station", varNames, value=TRUE)) == FALSE)]
            climVarList = list()
            for (k in 1:length(climVars)) {
                  length(climVarList) = k
                  vDescript = varList[[grep(paste("^", climVars[k], "$", sep=""), varNames)]]$getDescription()
                  vUnits = varList[[grep(paste("^", climVars[k], "$", sep=""), varNames)]]$getUnitsString()
                  climVarList[[k]] = list("Description" = vDescript, "Units" = vUnits)
            }
            names(climVarList) <- climVars
            var.list <- list("Station_Ids" = stationIds, "LonLatCoords" = lonLatList, "Altitudes" = altList, "Times" = timeList, "ClimVars" = climVarList)
      } else {
            # datasetType = "GridDataset"
            var.list <- list()
            for (i in 1:length(varNames)) {
                  dataVar <- gds$getDataVariable(varNames[i])
                  description <- dataVar$getDescription()
                  varName <- dataVar$getShortName()
                  dataType <- dataVar$getDataType()$toString()
                  dimensions <- unlist(strsplit(dataVar$getDimensionsString(), "\\s"))
                  units <- dataVar$getUnitsString()
                  grid <- gds$findGridByName(varName)
                  gridLevels <- as.numeric(unlist(strsplit(gsub("\\[|]|\\s","", grid$getLevels()$toString()), split=",")))
                  # Handles the non-existence of some levels for variables with vertical axis 
                  ## It is important to know the position of the dimension level (it cannot be guaranteed that it is the second)
                  levelDimIndex <- grep("^lev", dimensions)
                  if (length(gridLevels) > 0) {
                        v <- c()
                        for (k in 1:length(gridLevels) - 1) {
                              start <- rep(1, length(dimensions))
                              stride <- start
                              start[levelDimIndex] <- k
                              res <- tryCatch({ # Tries to read one value from each level
                                    grid$getVariable()$read(as.integer(start), as.integer(stride))$getStorage()
                              }, error = function(er) { # Handles the error when level does not exist
                                    err <- NA
                                    return(err)
                              })
                              v <- c(v, res)
                        }
                        noLevelInd <- which(is.na(v))
                  }     
                  gcs <- grid$getCoordinateSystem()
                  dim.list <- list()
                  for (j in 1:length(dimensions)) {
                        if(dimensions[j] == "member") {
                              axis <- gcs$getEnsembleAxis()
                        }
                        if(grepl("^run", dimensions[j])) {
                              axis <- gcs$getRunTimeAxis()
                        }
                        if (grepl("^time", dimensions[j])) {
                              axis <- gcs$getTimeAxis()
                        }
                        if(grepl("^lat", dimensions[j])) {
                              axis <- gcs$getLatAxis()
                        }
                        if(grepl("^lon", dimensions[j])) {
                              axis <- gcs$getLonAxis()
                        }
                        if(grepl("^lev|zeta", dimensions[j])) {
                              axis <- gcs$getVerticalAxis()
                        }
                        axisType <- axis$getAxisType()$toString()
                        # TODO: consider 2-D lon/lat dimension definition!!
                        if (grepl("^time|run|^lev", dimensions[j]) == FALSE) {
                              values <- axis$getCoordValues() 
                        }
                        if(grepl("^lev", dimensions[j])) {
                              if (length(noLevelInd) > 0) {
                                    values <- axis$getCoordValues()[-noLevelInd]
                              } else {
                                    values <- axis$getCoordValues()
                              }
                        }
                        if (grepl("^time", dimensions[j])) {
                              charDates <- unlist(strsplit(gsub("\\[|]|\\s", "", gcs$getTimes()$toString()), split = ","))
                              values <- strptime(charDates, format = "%Y-%m-%dT%H:%M:%SZ")
                              time.agg <- difftime(values[2], values[1], units = "hours")
                        }
                        if(grepl("^run", dimensions[j])) {
                              values <- strptime(gsub("\\[|]|\\s", "", unlist(strsplit(axis$getCalendarDates()$toString(), ","))), format = "%Y-%m-%dT%H:%M:%SZ")
                        }
                        dim.list[[j]] <- list("Type" = axisType, "Units" = axis$getUnitsString(), "Values" = values)
                  }
                  names(dim.list) <- dimensions
                  dimShape <- rep(NA, length(dim.list))
                  for (h in 1:length(dim.list)) {
                        dimShape[h] <- length(dim.list[[h]]$Values)
                  }
                  var.list[[i]] <- list("Description" = description, "DataType" = dataType, "Units" = units, "TimeStep" = time.agg, "Dimensions" = dim.list)
            }
            names(var.list) <- varNames
      }
      gds$close()
      return(var.list)
}
# End