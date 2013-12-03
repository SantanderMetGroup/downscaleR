dataInventory.ASCII <- function(dataset, rs) {
      dataset <- dataset
      rs <- rs
      lf <- list.files(dataset, full.names = TRUE)
      stations <- read.csv(lf[grep("stations", lf)], strip.white = TRUE)
      vars <- read.csv(lf[grep("variables", lf)], strip.white = TRUE, colClasses = "character")
      # VARS
      var.info <- list("varIDs" = vars[ ,1], "varNames" = vars[ ,2], "units" = vars[ ,3], "missing.code" = vars[ ,4])
      # DATES
      timeString <- read.csv(lf[grep(paste("^", vars[ ,1][1], "\\.*", sep = ""), list.files(dataset))], colClasses = "character")[ ,1]
      if (nchar(timeString[1]) == 8) {
            timeDates <- strptime(timeString, "%Y%m%d")  
      } else {
            timeDates <- strptime(timeString, "%Y%m%d%H")
      }
      rm(timeString)
      timeAxis <- list("startDate" = min(timeDates), "endDate" = max(timeDates), "timeStep" = difftime(timeDates[2], timeDates[1], units = "h"))    
      rm(timeDates)
      station.info <- list("stationIDs" = as.character(stations[ ,1]), "stationNames" = as.character(stations[ ,2]), "Altitude" = stations[ ,5], "LonLatCoords" = stations[ ,3:4], "timeAxis" = timeAxis)
      rm(timeAxis)
      info <- list("Stations" = station.info, "Variables" = var.info, "perc.missing" = NULL)
      if (isTRUE(rs)) {
            na.perc <- function(x) {round(length(which(is.na(x))) * 100 / length(x), 1)}
            aux.mat <- matrix(ncol = nrow(vars), nrow = length(station.info$stationIDs), dimnames = list(station.info$stationIDs, vars$ID))
            for (i in 1:nrow(vars)) {
                  var <- read.csv(lf[grep(paste("^", vars$ID[i], "\\.*", sep = ""), list.files(dataset))], na.strings = var.info$missing.code[i])[ ,-1]
                  aux.mat[ ,i] <- apply(var, 2, na.perc)
            }
            info <- list("Stations" = station.info, "Variables" = var.info, "perc.missing" = aux.mat)
      }
      return(info)
}
# End
# Example wiki [https://github.com/SantanderMetGroup/downscaling/wiki/Data-inventory]
# inventory.csv <- dataInventory(dataset = "datasets/observations/GSN_Iberia/", return.stats = TRUE)
# str(inventory.csv)
# inventory.csv$perc.missing
# which(inventory.csv$perc.missing == max(inventory.csv$perc.missing), arr.ind = TRUE)
# inventory.csv$Stations$stationNames[5]


