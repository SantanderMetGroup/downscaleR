#' Make data inventory from a station dataset in standard ASCII format
#' 
#' @param dataset path to the directory containng the dataset
#' @param rs 
#' @author J Bedia \email{joaquin.bedia@@gmail.com}
#' @keywords internal

dataInventory.ASCII <- function(dataset, rs) {
      lf <- list.files(dataset, full.names = TRUE)
      stations <- read.csv(lf[grep("stations", lf, ignore.case = TRUE)], strip.white = TRUE, stringsAsFactors = FALSE)
      vars <- read.csv(lf[grep("variables", lf, ignore.case = TRUE)], strip.white = TRUE, colClasses = "character")
      # Var info 
      var.info <- list("variable" = vars[ ,grep("variable", names(vars), ignore.case = TRUE)], "longname" = vars[ ,grep("longname", names(vars), ignore.case = TRUE)], "unit" = vars[ ,grep("unit", names(vars), ignore.case = TRUE)], "missing.code" = vars[ ,grep("missing_code", names(vars), ignore.case = TRUE)])
      var.info <- do.call("cbind.data.frame", var.info)
      # Station info
      timeString <- read.csv(lf[grep(paste("^", vars[ ,1][1], "\\.*", sep = ""), list.files(dataset))], colClasses = "character")[ ,1]
      if (nchar(timeString[1]) == 8) {
            timeDates <- strptime(timeString, "%Y%m%d")  
      } else {
            timeDates <- strptime(timeString, "%Y%m%d%H")
      }
      rm(timeString)
      timeAxis <- list("startDate" = min(timeDates), "endDate" = max(timeDates), "timeStep" = difftime(timeDates[2], timeDates[1], units = "h"))    
      rm(timeDates)
      station_id <- as.character(stations[ ,grep("station_id", names(stations), ignore.case = TRUE)])
      lon <- stations[ ,grep("^longitude$", names(stations), ignore.case = TRUE)]
      lat <- stations[ ,grep("^latitude$", names(stations), ignore.case = TRUE)]
      LonLatCoords <- cbind(lon, lat)
      rownames(LonLatCoords) <- station_id
      rm(lon, lat)
      other.metadata <- as.list(stations[ ,-pmatch(c("station_id", "longitude", "latitude"), names(stations))])
      station.info <- list("station_id" = station_id, "LonLatCoords" = LonLatCoords, "times" = timeAxis, "other.metadata" = other.metadata)
      rm(station_id, timeAxis, LonLatCoords, other.metadata)
      info <- list("Stations" = station.info, "Variables" = var.info, "Summary.stats" = NULL)
      if (isTRUE(rs)) {
            na.perc <- function(x) {round(length(which(is.na(x))) * 100 / length(x), 1)}
            aux.mat <- matrix(ncol = nrow(vars), nrow = length(station.info$station_id), dimnames = list(station.info$station_id, vars[ ,grep("variable", names(vars), ignore.case = TRUE)]))
            stats.list <- list("missing.percent" = aux.mat, "min" = aux.mat, "max" = aux.mat, "mean" = aux.mat)
            for (i in 1:nrow(vars)) {
                  var <- read.csv(lf[grep(paste("^", var.info[i,1], "\\.*", sep = ""), list.files(dataset))], na.strings = var.info$missing.code[i])[ ,-1]
                  stats.list[[1]][ ,i] <- apply(var, 2, na.perc)
                  stats.list[[2]][ ,i] <- apply(var, 2, min, na.rm = TRUE)
                  stats.list[[3]][ ,i] <- apply(var, 2, max, na.rm = TRUE)
                  stats.list[[4]][ ,i] <- apply(var, 2, mean, na.rm = TRUE)
                  rm(var)
            }
            info <- list("Stations" = station.info, "Variables" = var.info, "Summary.stats" = stats.list)
      }
      return(info)
}
# End

# # Example wiki [https://github.com/SantanderMetGroup/downscaling/wiki/Data-inventory]
# inventory.csv <- dataInventory(dataset = "datasets/observations/GSN_Iberia/", return.stats = TRUE)
# str(inventory.csv)
# inventory.csv$perc.missing
# which(inventory.csv$perc.missing == max(inventory.csv$perc.missing), arr.ind = TRUE)
# inventory.csv$Stations$stationNames[5]

